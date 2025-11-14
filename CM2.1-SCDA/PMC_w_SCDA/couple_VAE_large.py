import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F
from timm.layers import trunc_normal_, DropPath
from demo_convNext import Block as conv_block
from demo_convNext import LayerNorm

class CrossAttention(nn.Module):
    def __init__(self,in_chans_x,in_chans_y,hidden_chans,nH,dropout=0.):
        super().__init__()
        self.hidden_chans = hidden_chans
        self.nH = nH
        self.hd = hidden_chans // nH
        assert hidden_chans % nH == 0
       
        self.Q = nn.Linear(in_chans_x,hidden_chans,bias=False)
        self.K = nn.Linear(in_chans_y,hidden_chans,bias=False)
        self.V = nn.Linear(in_chans_y,hidden_chans,bias=False)
        self.softmax = nn.Softmax(dim=-1)
        self.dropout_p = dropout
        self.norm = nn.LayerNorm(hidden_chans)

    def forward(self,x,y = None):
        if (y == None): y = x

        B, Cx, H, W = x.shape
        _, Cy, _, _ = y.shape

        # 确保张量连续性，避免步长不匹配
        x = x.permute(0,2,3,1).contiguous().reshape(B,H*W,Cx)
        y = y.permute(0,2,3,1).contiguous().reshape(B,H*W,Cy)
        
        ## B,H*W,hidden_chans -- B,H*W,nH,hd -- B,nH,H*W,hd
        Qs = self.Q(x)

        Qh = Qs.view(B,H*W,self.nH,self.hd).transpose(1,2).contiguous()
        Kh = self.K(y).view(B,H*W,self.nH,self.hd).transpose(1,2).contiguous()
        Vh = self.V(y).view(B,H*W,self.nH,self.hd).transpose(1,2).contiguous()

        attn = F.scaled_dot_product_attention(Qh,Kh,Vh,attn_mask=None,
                  dropout_p= self.dropout_p, is_causal=False)
        
        # 确保注意力输出连续
        attn = attn.transpose(1,2).contiguous().reshape(B,H*W,self.hidden_chans)

        out = self.norm(attn + Qs)
        return out.reshape(B,H,W,self.hidden_chans).permute(0,3,1,2).contiguous()

class ConvNeXt(nn.Module):
    def __init__(self, in_chans=3, depths=[2,2,6,2], dims=[48,96,192,384],drop_path_rate=0.,
            layer_scale=1e-6):
        super().__init__()
        self.layers = len(depths)
        self.downsample_layers = nn.ModuleList()

        stem = nn.Sequential(nn.Conv2d(in_chans,dims[0],kernel_size=3,stride=1,padding=1),
                LayerNorm(dims[0]))
        self.downsample_layers.append(stem)

        for i in range(self.layers-1):
            self.downsample_layers.append(
               nn.Sequential(LayerNorm(dims[i]),
               nn.Conv2d(dims[i],dims[i+1],kernel_size=2,stride=2)))

        self.stages = nn.ModuleList()
        dp_rates = [x.item() for x in torch.linspace(0, drop_path_rate, sum(depths))]
        cur = 0
        for i in range(self.layers):
            self.stages.append(
              nn.Sequential( *[conv_block(dim=dims[i],drop_path=dp_rates[cur + j],
              layer_scale=layer_scale) for j in range(depths[i])]) )
            cur += depths[i]

        self.norm = LayerNorm(dims[-1],eps=1e-6)
        self.apply(self._init_weights)

    def _init_weights(self,m):
        if isinstance(m, (nn.Conv2d, nn.Linear)):
            trunc_normal_(m.weight,std=.02)
            nn.init.constant_(m.bias, 0)

    def forward(self,x):
        for i in range(self.layers):
            x = self.downsample_layers[i](x)
            x = self.stages[i](x)
        return self.norm(x)

class VAE_Encoder(nn.Module):
    def __init__(self,latent_dim=1024):
        super().__init__()
        self.cross_size = (16,32)
        self.atm_layer = ConvNeXt(in_chans=4,depths=[2,2,6],dims=[24,48,96],drop_path_rate=0.02) 
        self.ocn_layer = ConvNeXt(in_chans=4,depths=[2,2,6,2],dims=[24,48,96,192],drop_path_rate=0.05) 
        self.cross_atm = CrossAttention(in_chans_x=96,in_chans_y=192,hidden_chans=96,nH=8)
        self.cross_ocn = CrossAttention(in_chans_x=192,in_chans_y=96,hidden_chans=192,nH=8)
        self.conv_mixer = nn.Sequential(*[conv_block(dim=96+192,drop_path=0.,) for _ in range(2)], 
                           LayerNorm(96+192),
                           nn.Conv2d(96+192,192*2,kernel_size=3,stride=2,padding=1))

        self.out = nn.Sequential(LayerNorm(192*2), 
                      nn.GELU(),
                      nn.Conv2d(192*2,latent_dim,kernel_size=3,stride=1,padding=1))

        self.mlp = nn.Sequential(
                         nn.Flatten(2),
                         nn.Linear(8*16,64),
                         nn.GELU(),
                         nn.Linear(64,2))

    def forward(self,atm,ocn):
        atm_encode = self.atm_layer(atm)
        ocn_encode = self.ocn_layer(ocn)

        # 确保插值后的张量连续
        atm_encode = F.interpolate(atm_encode,size=self.cross_size,align_corners=True,mode='bilinear').contiguous()
        ocn_encode = F.interpolate(ocn_encode,size=self.cross_size,align_corners=True,mode='bilinear').contiguous()

        atm_attn = self.cross_atm(atm_encode,ocn_encode)
        ocn_attn = self.cross_ocn(ocn_encode,atm_encode)

        couple_mode = torch.cat([atm_attn,ocn_attn],dim=1)
        couple_mode = self.conv_mixer(couple_mode)
        latent = self.out(couple_mode)
        return  self.mlp(latent)

class Upsample(nn.Module):
    def __init__(self,in_chans,out_chans):
        super().__init__()
        self.up = nn.Sequential(
                nn.Upsample(scale_factor=2,mode='bilinear',align_corners=True),
                nn.Conv2d(in_chans,out_chans,kernel_size=3,padding=1))
    def forward(self,x):
        return self.up(x).contiguous()  # 确保输出连续

class VAE_Decoder(nn.Module):
    def __init__(self,latent_dim=1024):
        super().__init__()

        ##expand the spatial dim
        self.mlp = nn.Sequential(
                   nn.Linear(1,64),
                   nn.GELU(),
                   nn.Linear(64,16*8))

        # choose embed channel to control decoder size
        self.embed_ch = 1024   # 你可以试 512/1024/2048，看 trade-off

        # projection from latent_dim -> embed_ch
        self.pre_proj = nn.Conv2d(latent_dim, self.embed_ch, kernel_size=1)

        # conv_embed consumes embed_ch now
        self.conv_embed = Upsample(self.embed_ch, 192*4)

        self.cross_atm = CrossAttention(192*2,192*2,192,nH=8)
        self.cross_ocn = CrossAttention(192*2,192*2,192,nH=8)

        atm_modules = []
        in_chans, out_chans = 192,192
        for i in range(3):
            atm_modules.append(conv_block(in_chans)) 
            out_chans = in_chans // 2
            if i < (3-1):
                atm_modules.append(Upsample(in_chans,out_chans))
            in_chans = out_chans
        self.atm_decoder = nn.Sequential(*atm_modules)

        ocn_modules = []
        in_chans, out_chans = 192,192
        for i in range(4):
            ocn_modules.append(conv_block(in_chans)) 
            out_chans = in_chans // 2
            if i < (4-1):
                ocn_modules.append(Upsample(in_chans,out_chans))
            in_chans = out_chans
        self.ocn_decoder = nn.Sequential(*ocn_modules)
        ##atm_attn B,192,16,32 -- B,96,32,64 -- B,48,64,128 -- B,3,90,144
        ##ocn_attn B,192,16,32 -- B,96,32,64 -- B,48,64,128 -- B,24,128,256 -- B,3,200,360
        self.atm_interp =nn.Sequential(
            nn.Conv2d(48, 24, kernel_size=3, stride=1, padding=1),
            nn.GELU(),
            nn.Upsample(size=(90, 144), mode='bilinear', align_corners=True),
            nn.Conv2d(24, 4, kernel_size=3, stride=1, padding=1))

        self.ocn_interp =nn.Sequential(
            nn.Conv2d(24,12, kernel_size=3, stride=1, padding=1),
            nn.GELU(),
            nn.Upsample(size=(200, 360), mode='bilinear', align_corners=True),
            nn.Conv2d(12, 4, kernel_size=3, stride=1, padding=1))

    def forward(self,x):
        x = self.mlp(x)
        B,latent_dim,N = x.shape

        # 使用contiguous确保内存连续性
        x = x.contiguous().view(B,latent_dim,8,16)
        ###B,384*2,16,32
        # 假设 x 形状已是 (B, latent_dim, 8, 16)
        x = self.pre_proj(x)     # -> (B, embed_ch, 8, 16)
        x = self.conv_embed(x)
        _,C,_,_ = x.shape
        atm_dec = x[:,:C//2,:,:].contiguous()
        ocn_dec = x[:,C//2:,:,:].contiguous()

        atm_attn = self.cross_atm(atm_dec,ocn_dec)
        ocn_attn = self.cross_ocn(ocn_dec,atm_dec)

        atm_prm = self.atm_decoder(atm_attn)
        ocn_prm = self.ocn_decoder(ocn_attn)

        return self.atm_interp(atm_prm), self.ocn_interp(ocn_prm)


# 推荐保存在模型类或模块顶端
LOGVAR_MIN = -30.0
LOGVAR_MAX = 6.0   # 推荐从 10 开始，若仍不稳再降到 6
MEAN_ABS_CLAMP = 20.0  # 可选：限制 mean 的绝对值
MEAN_L2_REG = 1e-6     # 可选：对 mean 的 L2 正则强度（0 为不开）
GRAD_CLIP_DEFAULT = 2.0

class VAE(nn.Module):
    def __init__(self,latent_dim=2048):
        super().__init__()
        self.encoder = VAE_Encoder(latent_dim=latent_dim)
        self.decoder = VAE_Decoder(latent_dim=latent_dim)

    def safe_reparam(self, mean, logvar):
        # clamp logvar BEFORE using it to compute std
        logvar_clamped = torch.clamp(logvar, min=LOGVAR_MIN, max=LOGVAR_MAX)
        std = torch.exp(logvar_clamped * 0.5)
        eps = torch.randn_like(std)
        z = eps * std + mean
        return z, mean, logvar_clamped

    def encode(self,x_atm,x_ocn):
        encoded = self.encoder(x_atm,x_ocn)
        mean   = encoded[:,:,0]
        logvar = encoded[:,:,1]

        ##reparameterization
        if self.training:
            z, mean, logvar = self.safe_reparam(mean, logvar)
            #z = mean
        else:
            # 评估时直接使用 clamp 后的 mean/logvar（避免 eval 时异常）
            logvar = torch.clamp(logvar, min=LOGVAR_MIN, max=LOGVAR_MAX)
            z = mean
        return [z,mean,logvar]

    def decode(self,z):
        ##z shape is: (B,latent_dim)
        z = z.unsqueeze(-1)
        return self.decoder(z) 

    def forward(self,atm,ocn):
        z,mean,logvar = self.encode(atm,ocn)
        decoded = self.decode(z)
        return decoded, mean, logvar