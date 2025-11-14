import torch
import torch.nn as nn
import torch.nn.functional as F
from timm.layers import trunc_normal_, DropPath

class LayerNorm(nn.Module):
    def __init__(self,dim,eps=1e-6,data_format="C_First"):
        super().__init__()
        self.weight = nn.Parameter(torch.ones(dim))
        self.bias = nn.Parameter(torch.ones(dim))
        self.eps = eps
        self.data_format = data_format
        self.dim = (dim,)
    def forward(self,x):
        if self.data_format == "C_Last":
            ##[N,H,W,C]
            return F.layer_norm(x,self.dim,self.weight,self.bias,self.eps)
        elif self.data_format == "C_First":
            ##[N,C,H,W]
            u = x.mean(1,keepdim=True)
            s = (x - u).pow(2).mean(1, keepdim=True)
            x = (x - u) / torch.sqrt(s + self.eps)
            assert x.shape[1] == self.weight.shape[0], f"Expected {self.weight.shape[0]} channels, got {x.shape[1]}"
            x = self.weight[:,None,None]*x  + self.bias[:,None,None]
            return x

class Block(nn.Module):
    def __init__(self, dim, drop_path=0., layer_scale=1e-6):
        super().__init__()
        ##depth-wise conv
        self.dwconv = nn.Conv2d(dim, dim, kernel_size=7, padding=3, groups=dim)
        self.norm = LayerNorm(dim,data_format="C_Last")
        ##point-wise conv
        self.pwconv1 = nn.Linear(dim, 4*dim)
        self.act = nn.GELU()
        self.pwconv2 = nn.Linear(4*dim,dim)
        self.gamma = nn.Parameter(layer_scale * torch.ones((dim)), 
           requires_grad=True) if layer_scale >0 else None
        self.drop_path = DropPath(drop_path) if drop_path > 0 else nn.Identity()

    def forward(self,x):
        shortcut = x
        x = self.dwconv(x)
        x = x.permute(0,2,3,1) ##(N,C,H,W) -- (N,H,W,C)
        x = self.norm(x)
        x = self.pwconv1(x)
        x = self.act(x)
        x = self.pwconv2(x)
        if self.gamma is not None:
            x = self.gamma * x
        x = x.permute(0,3,1,2) ##N,C,H,W
        x = shortcut + self.drop_path(x)
        return self.drop_path(x) + shortcut
 
class ConvNeXt(nn.Module):
    def __init__(self, in_chans=3, depths=[3,3,9,3], dims=[96,192,384,768],drop_path_rate=0.,
            layer_scale=1e-6):
        super().__init__()
        self.layers = len(depths)
        self.downsample_layers = nn.ModuleList()

        stem = nn.Sequential(nn.Conv2d(in_chans,dims[0],kernel_size=4,stride=4),
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
              nn.Sequential( *[Block(dim=dims[i],drop_path=dp_rates[cur + j],
              layer_scale=layer_scale) for j in range(depths[i])]) )
            cur += depths[i]

        self.norm = nn.LayerNorm(dims[-1],eps=1e-6)
        self.apply(self._init_weights)

    def _init_weights(self,m):
        if isinstance(m, (nn.Conv2d, nn.Linear)):
            trunc_normal_(m.weight,std=.02)
            nn.init.constant_(m.bias, 0)

    def forward(self,x):
        for i in range(self.layers):
            x = self.downsample_layers[i](x)
            x = self.stages[i](x)
        return self.norm(x.mean([-2,-1])) ##gloabl average pooling


if __name__ == "__main__":
    a = ConvNeXt()
    a = a.to('cuda')
    from torchsummary import summary
    summary(a,input_size=(3,200,360),device='cuda')
    ##the output shape is batch_size,768


