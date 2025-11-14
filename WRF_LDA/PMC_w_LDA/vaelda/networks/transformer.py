import torch.nn as nn
import torch
import copy
from functools import partial
from einops import rearrange
from timm.models.layers import DropPath, to_2tuple, trunc_normal_
from fairscale.nn.checkpoint.checkpoint_activations import checkpoint_wrapper
from torchvision import utils as vutils
from vaelda.networks.utils.swinblock import SwinTransformerBlock
# from torch.distributed.algorithms._checkpoint._checkpoint_wrapper import checkpoint_wrapper


class PatchEmbed(nn.Module):
    r""" Image to Patch Embedding【切块】图像 → Token
    Args:
        img_size (int): 输入图像大小.  Default: 224.
        patch_size (int): Patch token size. Default: 4.
        in_chans (int): 输入图像通道数. Default: 3.
        embed_dim (int): 输出维度. Default: 96.
        norm_layer (nn.Module, optional): Normalization layer. Default: None
    """

    def __init__(self, img_size=[256, 256], patch_size=4, stride=[4,4],
                 in_chans=4, embed_dim=96, norm_layer=None):
        super().__init__()
        patches_resolution = [img_size[0] // stride[0], img_size[1] // stride[1]]
        # 计算输出的 patch 数量（纵向 × 横向），用于后续形状计算。
        self.img_size = img_size
        self.patch_size = patch_size
        self.patches_resolution = patches_resolution
        self.num_patches = patches_resolution[0] * patches_resolution[1]

        # 记录图像信息，num_patches 是总的 patch 数量。
        self.in_chans = in_chans
        self.embed_dim = embed_dim

        # 保存输入输出通道数。使用卷积代替显式的 patch 切分与线性投影操作
        self.proj = nn.Conv2d(in_chans, embed_dim, kernel_size=patch_size, stride=stride)
        # 将每个 patch 映射为一个 embed_dim 维的 token。如果设置了归一化层，则添加到网络中。
        if norm_layer is not None:
            self.norm = norm_layer(embed_dim)
        else:
            self.norm = None

    def forward(self, x):
        B, C, H, W = x.shape
        # FIXME look at relaxing size constraints
        # 提取 batch 大小、通道数、高度、宽度。
        assert H == self.img_size[0] and W == self.img_size[1], \
            f"Input image size ({H}*{W}) doesn't match model ({self.img_size[0]}*{self.img_size[1]})."
        # 检查输入尺寸是否符合预设值。
        # self.proj(x)：使用卷积将图像映射为 patch token。
        # .flatten(2)：将 H×W 合并成一个维度（即 Ph×Pw）。
        # .transpose(1, 2)：调整形状为 [B, num_patches, embed_dim]，便于输入 Transformer。
        x = self.proj(x).flatten(2).transpose(1, 2)  # B Ph*Pw C
        if self.norm is not None:
            x = self.norm(x)
        return x                  # 输出 [B, num_patches, embed_dim]

    # flops 方法
    # 计算卷积操作的 FLOPs，供性能分析用。
    def flops(self):
        Ho, Wo = self.patches_resolution
        flops = Ho * Wo * self.embed_dim * self.in_chans * (self.patch_size[0] * self.patch_size[1])
        # 归一化层的 FLOPs。
        if self.norm is not None:
            flops += Ho * Wo * self.embed_dim
        return flops

class PatchMerging(nn.Module):
    r""" Patch Merging Layer【下采样】空间降采样（通道加倍）.
    Args:
        input_resolution (tuple[int]): 输入特征的高宽.
        dim (int): 每个 token 的通道数.
        norm_layer (nn.Module, optional): Normalization layer.  Default: nn.LayerNorm
    """
    # input_resolution: 输入特征的高宽。
    # dim: 每个 token 的通道数。
    # norm_layer: 使用的归一化方法，默认是 LayerNorm。
    def __init__(self, input_resolution, dim, norm_layer=nn.LayerNorm):
        super().__init__()
        # 记录输入的维度信息。
        self.input_resolution = input_resolution
        self.dim = dim
        # 对 2x2 patch 合并后得到的 4×C 输入进行线性映射，输出维度为 2×C。
        self.reduction = nn.Linear(4 * dim, 2 * dim, bias=False)
        self.norm = norm_layer(4 * dim)

    def forward(self, x):
        """
        x: B, H*W, C
        """
        # # 输入形状是 [B, H, W, C]。
        H, W = self.input_resolution
        B, H, W, C = x.shape
        # 必须能整除 2，否则无法进行 2x2 patch 合并。
        assert H % 2 == 0 and W % 2 == 0, f"x size ({H}*{W}) are not even."

        # x = x.view(B, H, W, C)

        # 取出每 2x2 网格中的四个角 patch，形状都是 [B, H/2, W/2, C]。
        x0 = x[:, 0::2, 0::2, :]  # B H/2 W/2 C
        x1 = x[:, 1::2, 0::2, :]  # B H/2 W/2 C
        x2 = x[:, 0::2, 1::2, :]  # B H/2 W/2 C
        x3 = x[:, 1::2, 1::2, :]  # B H/2 W/2 C
        # 拼接四个 patch 的通道维，结果为 [B, H/2, W/2, 4*C]。
        x = torch.cat([x0, x1, x2, x3], -1)  # B H/2 W/2 4*C
        # reshape 为 [B, num_new_patches, 4*C]。
        x = x.view(B, -1, 4 * C)  # B H/2*W/2 4*C

        # 对每个 patch 的 4×C 维向量进行归一化。
        x = self.norm(x).view(B, H//2, W//2, 4*C)
        # 通过线性层将其降维到 [B, H/2, W/2, 2*C]。
        x = self.reduction(x)

        return x


class PatchExpand(nn.Module):
    """
    上采样模块 上采样恢复图像空间
    dim: 输入通道数
    dim_scale: 扩大维度的倍数（默认是2，意味着2倍上采样）
    norm_layer: 归一化层
    """
    def __init__(self, dim, dim_scale=2, norm_layer=nn.LayerNorm):
        super().__init__()
        # 如果需要上采样，则通过线性层将 dim → 2×dim。
        # 否则，直接使用 Identity（无变化）。
        self.dim = dim
        self.dim_scale = dim_scale
        self.expand = nn.Linear(dim, 2*dim, bias=False) if dim_scale==2 else nn.Identity()
        # 上采样后维度会变成 dim//dim_scale，因此归一化的输入是下采样后的通道数。
        self.norm = norm_layer(dim // dim_scale)

    def forward(self, x):
        """
        x: B, H*W, C
        """
        # 输入是 [B, H*W, C]，这是 Transformer 输出的格式。
        B, L, C = x.shape

        # 根据 dim 推断原图像 H 和 W
        H = W = int(L ** 0.5)
        assert H * W == L, f"PatchExpand 输入不是正方形，x.shape={x.shape}"

        # 线性变换维度扩大。
        x = self.expand(x)  # # B, L, 2C
        # 这行代码实际上有 bug，因为 x 此时是 [B, N, C]，需要 reshape 成 [B, H, W, C] 之前，
        # 你应该先知道 H 和 W。
        # 但是此处你代码未做 reshape，意味着你原本使用 x = x.view(B, H, W, C) 是注释掉的。
        # 后面的 rearrange 会假设当前就是 [B, H, W, C]。
        # B, H, W, C = x.shape

        x = x.view(B, H, W, -1)  # # B, H, W, 2C
        # 这是关键：将通道维度展开为 spatial 上的尺寸，也就是：
        # 原通道数 C = 4 × c，被 reshape 成 2 × 2 × c
        # 然后 h → h×2，w → w×2，即空间维度放大 2 倍，通道变为 C//4
        # 所以这是空间上的上采样
        x = rearrange(x, 'b h w (p1 p2 c)-> b (h p1) (w p2) c', p1=2, p2=2)
        # x = x.view(B,-1,C//4)
        # 归一化后输出。输出维度：[B, H×2, W×2, C//4]
        x= self.norm(x)

        return x

class FinalPatchExpand_X4(nn.Module):
    """
    4倍上采样模块 上采样恢复图像空间
    input_resolution: 当前特征的高和宽
    dim: 当前通道数
    dim_scale: 上采样倍数，默认是 4
    """
    def __init__(self, input_resolution, dim, dim_scale=4, norm_layer=nn.LayerNorm):
        # 线性层将 dim → 16×dim，为后续 reshape 成 4×4 patch 做准备
        super().__init__()
        self.input_resolution = input_resolution
        self.dim = dim
        self.dim_scale = dim_scale
        self.expand = nn.Linear(dim, 16*dim, bias=False)
        # 输出通道仍设为 dim，加 norm。
        self.output_dim = dim
        self.norm = norm_layer(self.output_dim)

    def forward(self, x):
        """
        x: B, H*W, C
        """
        # 将每个 token 扩展为 16×dim 长度。
        H, W = self.input_resolution
        x = self.expand(x)
        # 确认 token 数和期望一致。
        B, L, C = x.shape
        assert L == H * W, "input feature has wrong size"

        # reshape 成图像空间的格式。
        x = x.view(B, H, W, C)
        # 将 C = 16×dim → (p1=4, p2=4, c)，即变成 4×4×dim 的 patch。
        # spatial 维度扩大为 H×4, W×4。
        x = rearrange(x, 'b h w (p1 p2 c)-> b (h p1) (w p2) c', p1=self.dim_scale, p2=self.dim_scale, c=C//(self.dim_scale**2))
        # 输出恢复为 [B, H×W, C] 格式，适合输出或可视化。
        x = x.view(B,-1,self.output_dim)
        x= self.norm(x)

        return x



class BasicLayer(nn.Module):
    """ A basic Swin Transformer layer for one stage【编码模块】逐层编码.
    Args:
        dim (int): Number of input channels.
        input_resolution (tuple[int]): Input resolution.
        depth (int): Number of blocks.
        num_heads (int): Number of attention heads.
        window_size (int): Local window size.
        mlp_ratio (float): Ratio of mlp hidden dim to embedding dim.
        qkv_bias (bool, optional): If True, add a learnable bias to query, key, value. Default: True
        qk_scale (float | None, optional): Override default qk scale of head_dim ** -0.5 if set.
        drop (float, optional): Dropout rate. Default: 0.0
        attn_drop (float, optional): Attention dropout rate. Default: 0.0
        drop_path (float | tuple[float], optional): Stochastic depth rate. Default: 0.0
        norm_layer (nn.Module, optional): Normalization layer. Default: nn.LayerNorm
        downsample (nn.Module | None, optional): Downsample layer at the end of the layer. Default: None
        use_checkpoint (bool): Whether to use checkpointing to save memory. Default: False.
    """

    def __init__(self, dim, input_resolution, depth, num_heads, window_size,
                 mlp_ratio=4.0, qkv_bias=True, qk_scale=None, drop=0.0, attn_drop=0.0,
                 drop_path=0.1, norm_layer=nn.LayerNorm, downsample=None,
                 use_checkpoint=False, pre_norm=True):

        super().__init__()
        self.dim = dim
        self.input_resolution = input_resolution
        self.depth = depth
        self.use_checkpoint = use_checkpoint

        # build blocks
        # 构建多个 Swin Block，交替设置 shift（实现 SW-MSA 和 SW-MSA with shift）。
        self.blocks = nn.ModuleList()
        for i in range(depth):
            block = SwinTransformerBlock(
                dim=dim,
                input_resolution=input_resolution,
                num_heads=num_heads,
                window_size=window_size,
                shift_size=0 if i%2==0 else window_size//2,
                mlp_ratio=4.,
                qkv_bias=True,
                qk_scale=None,
                drop=0.,
                attn_drop=0.,
                drop_path=0.,
            )

            # 支持 checkpoint 以节省内存
            if use_checkpoint:
                block = checkpoint_wrapper(block, offload_to_cpu=True)
            self.blocks.append(block)


        if downsample is not None:
            self.downsample = downsample(input_resolution, dim=dim//2, norm_layer=norm_layer)
        else:
            self.downsample = None

    def forward(self, x):
        # 先下采样（一般是 stage 的入口）
        if self.downsample is not None:
            x = self.downsample(x)
        # 依次通过 Swin Transformer Block。
        for blk in self.blocks:
            # if self.use_checkpoint:
            #     x = checkpoint.checkpoint(blk, x)
            # else:
            #     x = blk(x)
            x = blk(x)
        return x


class BasicLayer_up(nn.Module):
    """
        一个上采样模块：多个 SwinTransformerBlock + PatchExpand

        Args:
            dim: 输入特征维度
            input_resolution: 输入的空间尺寸，如 (H, W)
            depth: 本层包含多少个 SwinTransformerBlock
            upsample: 是否需要 PatchExpand 上采样（一般除了最后一层都需要）
    """

    def __init__(self, dim, input_resolution, depth, num_heads, window_size,
                 mlp_ratio=4.0, qkv_bias=True, qk_scale=None, drop=0.0, attn_drop=0.0,
                 drop_path=0.1, norm_layer=nn.LayerNorm, upsample=None,
                 use_checkpoint=False, pre_norm=True):

        super().__init__()
        self.dim = dim
        self.depth = depth
        self.use_checkpoint = use_checkpoint

        # 构建 SwinTransformerBlock
        self.blocks = nn.ModuleList([
            SwinTransformerBlock(
                dim=dim,
                input_resolution=input_resolution,
                num_heads=num_heads,
                window_size=window_size,
                shift_size=0 if (i % 2 == 0) else window_size // 2,  # 交替使用 shift
                mlp_ratio=mlp_ratio,
                qkv_bias=qkv_bias, qk_scale=qk_scale,
                drop=drop, attn_drop=attn_drop,
                drop_path=drop_path[i] if isinstance(drop_path, list) else drop_path,
                norm_layer=norm_layer,
            ) for i in range(depth)
        ])

        # 如果设置了 upsample（即不是最后一层），则加上 PatchExpand
        if upsample is not None:
            self.upsample = upsample(dim=dim, norm_layer=norm_layer)
        else:
            self.upsample = None

    def forward(self, x):
        # Transformer Block 前置（解码时先 refine 后恢复分辨率）
        """
        x: [B, H, W, C]
        返回：
            如果设置 upsample，则输出上采样后的结果；
            否则返回原尺寸的处理结果。
        """
        for blk in self.blocks:
            # if self.use_checkpoint:
            #     x = checkpoint.checkpoint(blk, x)
            # else:
            #     x = blk(x)
            x = blk(x)  # SwinTransformerBlock 中默认接受 [B, H, W, C]
        # 最后上采样
        if self.upsample is not None:
            # Flatten 前两维为 L，变形为 [B, L, C] 传给 PatchExpand
            B, H, W, C = x.shape
            x = x.view(B, H * W, C)
            x = self.upsample(x)
        return x


class Transformer_Encoder(nn.Module):
    """
    Encoder（Swin-like）模块
    接收：图像张量 -> 输出：编码器最深层特征
    """
    def __init__(self, img_size=(256, 256), patch_size=4,  stride=(4,4), in_chans=4,
                 embed_dim=96, depths=[2, 2, 6], num_heads=[3, 6, 12],
                 window_size=8, mlp_ratio=4., qkv_bias=True, qk_scale=None,
                 drop_rate=0., attn_drop_rate=0., drop_path_rate=0.1,
                 norm_layer=nn.LayerNorm, ape=False, patch_norm=True,
                 use_checkpoint=False, pre_norm=True, **kwargs):
        super().__init__()


        self.num_layers = len(depths)
        self.embed_dim = embed_dim
        self.ape = ape
        self.patch_norm = patch_norm
        self.num_features = int(embed_dim * 2 ** (self.num_layers - 1))
        self.num_features_up = int(embed_dim * 2)
        self.mlp_ratio = mlp_ratio

        # Patch embedding: split image into non-overlapping patches
        self.patch_embed = PatchEmbed(
            img_size=img_size, patch_size=patch_size, stride=stride, in_chans=in_chans, embed_dim=embed_dim,
            norm_layer=norm_layer if self.patch_norm else None)
        num_patches = self.patch_embed.num_patches
        patches_resolution = self.patch_embed.patches_resolution
        self.patches_resolution = patches_resolution

        # Optional absolute positional embedding
        if self.ape:
            self.absolute_pos_embed = nn.Parameter(torch.zeros(1, num_patches, embed_dim))
            trunc_normal_(self.absolute_pos_embed, std=.02)

        self.pos_drop = nn.Dropout(p=drop_rate)

        # stochastic depth decay
        dpr = [x.item() for x in torch.linspace(0, drop_path_rate, sum(depths))]

        # build encoder layers
        self.layers = nn.ModuleList()
        for i_layer in range(self.num_layers):
            layer = BasicLayer(dim=int(embed_dim * 2 ** i_layer),
                               input_resolution=(patches_resolution[0] // (2 ** i_layer),
                                                 patches_resolution[1] // (2 ** i_layer)),
                               depth=depths[i_layer],
                               num_heads=num_heads[i_layer],
                               window_size=window_size,
                               mlp_ratio=self.mlp_ratio,
                               qkv_bias=qkv_bias, qk_scale=qk_scale,
                               drop=drop_rate, attn_drop=attn_drop_rate,
                               drop_path=dpr[sum(depths[:i_layer]):sum(depths[:i_layer + 1])],
                               norm_layer=norm_layer,
                               downsample=PatchMerging if (i_layer > 0) else None,
                               use_checkpoint=use_checkpoint,
                               pre_norm=pre_norm)
            self.layers.append(layer)

        # final norm on deepest feature
        self.norm = norm_layer(self.num_features)
        self.apply(self._init_weights)

    def _init_weights(self, m):
        if isinstance(m, nn.Linear):
            trunc_normal_(m.weight, std=.02)
            if isinstance(m, nn.Linear) and m.bias is not None:
                nn.init.constant_(m.bias, 0)
        elif isinstance(m, nn.LayerNorm):
            nn.init.constant_(m.bias, 0)
            nn.init.constant_(m.weight, 1.0)

    def forward(self, x):
        """
        输入: x -> [B, C, H, W]
        输出: 最深层编码特征 x -> [B, H', W', C']
        不再返回中间层的 skip features
        """
        B = x.shape[0]
        x = self.patch_embed(x)     # B, L, C
        if self.ape:
            x = x + self.absolute_pos_embed
        x = self.pos_drop(x)

        # 转为图像空间 [B, H, W, C]
        x = x.view(B, self.patches_resolution[0], self.patches_resolution[1], -1)
        for i, layer in enumerate(self.layers):
            x = layer(x)

        x = self.norm(x)  # B H W C
        return x


class Transformer_Decoder(nn.Module):
    """
    Decoder（对称于 Encoder），不使用 skip connection：
    仅接收 bottleneck 输出并逐层上采样恢复空间分辨率。
    保留原来层级与上采样模块设计，但去除与 encoder 的拼接逻辑。
    """
    def __init__(self, img_size=(256, 256), patch_size=4, stride=(4,4), in_chans=4, num_classes=4,
                 embed_dim=96, depths=[2, 2, 2], depths_decoder=[1, 2], num_heads=[12, 6, 3],
                 window_size=8, mlp_ratio=4., qkv_bias=True, qk_scale=None,
                 drop_rate=0., attn_drop_rate=0., drop_path_rate=0.1,
                 norm_layer=nn.LayerNorm, ape=False, patch_norm=True,
                 use_checkpoint=False, pre_norm=True, final_upsample="expand_first", **kwargs):
        super().__init__()

        self.num_classes = num_classes
        self.num_layers = len(depths)
        self.embed_dim = embed_dim
        self.ape = ape
        self.patch_norm = patch_norm
        self.num_features = int(embed_dim * 2 ** (self.num_layers - 1))
        self.num_features_up = int(embed_dim * 2)
        self.mlp_ratio = mlp_ratio
        self.final_upsample = final_upsample

        # stochastic depth
        dpr = [x.item() for x in torch.linspace(0, drop_path_rate, sum(depths))]

        # 构建解码器层（上采样层序列）
        self.layers_up = nn.ModuleList()
        for i_layer in range(self.num_layers):
            # idx corresponds to encoder stage index if mirrored, but we do not use skip features
            layer_up = BasicLayer_up(dim=int(embed_dim * 2 ** (self.num_layers-1-i_layer)),
                            depth=depths[(self.num_layers-1-i_layer)],
                            num_heads=num_heads[(self.num_layers-1-i_layer)],
                            window_size=window_size,
                            input_resolution=[img_size[0]//stride[0]//2**(self.num_layers-1-i_layer), img_size[1]//stride[1]//2**(self.num_layers-1-i_layer)],
                            mlp_ratio=self.mlp_ratio,
                            qkv_bias=qkv_bias, qk_scale=qk_scale,
                            drop=drop_rate, attn_drop=attn_drop_rate,
                            drop_path=dpr[sum(depths[:(self.num_layers-1-i_layer)]):sum(depths[:(self.num_layers-1-i_layer) + 1])],
                            norm_layer=norm_layer,
                            upsample=PatchExpand if (i_layer < self.num_layers - 1) else None,
                            use_checkpoint=use_checkpoint,
                            pre_norm=pre_norm)

            self.layers_up.append(layer_up)

        self.norm_up= norm_layer(self.embed_dim)
        self.apply(self._init_weights)

    def _init_weights(self, m):
        if isinstance(m, nn.Linear):
            trunc_normal_(m.weight, std=.02)
            if isinstance(m, nn.Linear) and m.bias is not None:
                nn.init.constant_(m.bias, 0)
        elif isinstance(m, nn.LayerNorm):
            nn.init.constant_(m.bias, 0)
            nn.init.constant_(m.weight, 1.0)

    def forward(self, x):
        """
        输入: bottleneck 特征 x -> [B, H', W', C]
        处理: 逐层上采样/变换
        输出: [B, L, C] 的 token（与原逻辑一致的输出格式）
        """
        for layer_up in self.layers_up:
            x = layer_up(x)

        x = self.norm_up(x)  # B H W C
        return x


class Layer(nn.Module):
    """
    基础 Transformer 层封装（共享结构）
    由若干 SwinTransformerBlock 叠加而成
    """
    def __init__(self, dim, depth, window_size, input_resolution,
                num_heads, mlp_ratio=4., qkv_bias=True,
                drop=0., attn_drop=0., drop_path=0.,
                norm_layer=nn.LayerNorm,
                use_checkpoint=False, pre_norm=True) -> None:
        super().__init__()
        self.dim = dim
        self.depth = depth
        self.save_attn = False

        self.blocks = nn.ModuleList()
        for i in range(depth):
            block = SwinTransformerBlock(
                dim=dim,
                input_resolution=input_resolution,
                num_heads=num_heads,
                window_size=window_size,
                shift_size=0 if i%2==0 else window_size//2,
                mlp_ratio=4.,
                qkv_bias=True,
                qk_scale=None,
                drop=drop,
                attn_drop=attn_drop,
                drop_path=drop_path[i] if isinstance(drop_path, list) else drop_path,
            )

            if use_checkpoint:
                block = checkpoint_wrapper(block, offload_to_cpu=True)
            self.blocks.append(block)

    def forward(self, x):
        for i, blk in enumerate(self.blocks):
            x = blk(x)
        return x


class Enc_net(nn.Module):
    """
    局部编码器网络：将输入一次性编码为高层特征。
    """
    def __init__(self, img_size=[256, 256], patch_size=4, stride=[4, 4], inchans_list=[4],
                 embed_dim=96, lgnet_embed=768, depths=[2,2,6], num_heads=[3,6,12],
                 window_size=8, inp_length=1, mlp_ratio=4, qkv_bias=True, drop_rate=0.,
                 attn_drop_rate=0., drop_path_rate=0., norm_layer=partial(nn.LayerNorm, eps=1e-6),
                 patch_norm=False, use_checkpoint=False, pre_norm=True):
        super().__init__()

        in_chans = sum(inchans_list)
        self.encoder = Transformer_Encoder(
            img_size=img_size,
            patch_size=patch_size,
            stride=stride,
            in_chans=in_chans,
            embed_dim=embed_dim,
            depths=depths,
            num_heads=num_heads,
            window_size=window_size,
            inp_length=inp_length,
            mlp_ratio=mlp_ratio,
            qkv_bias=qkv_bias,
            drop_rate=drop_rate,
            attn_drop_rate=attn_drop_rate,
            drop_path_rate=drop_path_rate,
            norm_layer=norm_layer,
            ape=True,
            patch_norm=patch_norm,
            use_checkpoint=use_checkpoint,
            pre_norm=pre_norm
        )

        # 将 encoder 最深特征线性投影到 lgnet 所需维度
        self.proj = nn.Linear(embed_dim * (2 ** (len(depths) - 1)), lgnet_embed)

    def forward(self, x):
        # x: (B, C, H, W)
        x = self.encoder(x)  # shape: (B, H', W', C')
        assert isinstance(x, torch.Tensor)
        B, H, W, C = x.shape

        x = x.view(B, -1, C)  # Flatten spatial
        x = self.proj(x)      # Linear projection to lgnet_embed
        return x


class Dec_net(nn.Module):
    """
    多头解码器网络（不再使用 skip connections）。
    将全局特征分配到各个分支并解码为图像形式。
    """
    def __init__(self, img_size=[256, 256], patch_size=4, stride=[4,4], outchans_list=[4], channel_num=4, embed_dim=96, lgnet_embed=768, depths=[2, 2, 2],
                 num_heads=[12, 6, 3], window_size=8, mlp_ratio=4, qkv_bias=True, drop_rate=0.,
                 attn_drop_rate=0., drop_path_rate=0., norm_layer=partial(nn.LayerNorm, eps=1e-6),
                 patch_norm=False,use_checkpoint=False, use_mlp=False, pre_norm=True):
        super().__init__()
        self.img_size = img_size
        self.patch_size = patch_size
        self.stride = stride
        self.embed_dim = embed_dim
        self.depths = depths
        self.outchans_list = outchans_list

        self.dec_list = nn.ModuleList()
        self.final_proj_list = nn.ModuleList()
        # 每个输出分支对应一个 Transformer_Decoder 与随后反卷积层
        for i in range(len(self.outchans_list)):
            self.dec_list.append(Transformer_Decoder(img_size=img_size, patch_size=patch_size, stride=stride, in_chans=self.outchans_list[i],
                                                     embed_dim=embed_dim, depths=depths, num_heads=num_heads, window_size=window_size,
                                                     mlp_ratio=mlp_ratio, qkv_bias=qkv_bias, drop_rate=drop_rate, attn_drop_rate=attn_drop_rate,
                                                     drop_path_rate=drop_path_rate, norm_layer=norm_layer, ape=False, patch_norm=patch_norm,
                                                     use_checkpoint=use_checkpoint, pre_norm=pre_norm))
            self.final_proj_list.append(nn.ConvTranspose2d(in_channels=embed_dim, out_channels=self.outchans_list[i],
                                                           kernel_size=patch_size, stride=stride))

        # 将 lgnet 输出投影回 decoder 所需的初始通道与分支数量
        self.proj = nn.Linear(lgnet_embed, embed_dim*2**(len(self.depths)-1)*len(self.outchans_list))


    def forward(self, x):
        # x: [B, L, lgnet_embed]
        data_proj = self.proj(x)
        # 拆分给各个分支
        data_split = torch.split(data_proj, self.embed_dim*2**(len(self.depths)-1), dim=-1)
        out_data_mean_list = []
        out_data_std_list = []
        for i in range(len(self.outchans_list)):
            # 每个解码器不再需要 downsample_list/skip connections
            out_data = self.dec_list[i](data_split[i]).permute(0,3,1,2)
            out_data = self.final_proj_list[i](out_data)
            out_data_mean_list.append(out_data[:, :out_data.shape[1]//2])
            out_data_std_list.append(out_data[:, out_data.shape[1]//2:])

        out_data_mean = torch.cat(out_data_mean_list, dim=1)
        out_data_std = torch.cat(out_data_std_list, dim=1)

        out_data = torch.cat((out_data_mean, out_data_std), dim=1)
        return out_data


class LG_net(nn.Module):
    """
    全局（LG_net）特征提取
    """
    def __init__(self, img_size=[256, 256], patch_size=4, in_chans=4, out_chans=4,
                 embed_dim=768, depths=(2, 2), num_heads=(6, 12),
                 window_size=8, mlp_ratio=4., qkv_bias=True,
                 drop_rate=0., attn_drop_rate=0., drop_path_rate=0.1,
                 norm_layer=partial(nn.LayerNorm, eps=1e-6), patch_norm=False,
                 use_checkpoint=False, pre_norm=True):
        super().__init__()

        self.num_layers = len(depths)
        self.embed_dim = embed_dim
        self.patch_norm = patch_norm
        # stage4输出特征矩阵的channels
        self.mlp_ratio = mlp_ratio
        self.window_size = window_size
        self.img_size = img_size
        self.patch_size = patch_size

        self.pos_drop = nn.Dropout(p=drop_rate)

        # stochastic depth
        dpr = [x.item() for x in torch.linspace(0, drop_path_rate, sum(depths))]  # stochastic depth decay rule

        # build layers
        self.layers = nn.ModuleList()
        for i_layer in range(self.num_layers):
            # 注意这里构建的stage和论文图中有些差异
            # 这里的stage不包含该stage的patch_merging层，包含的是下个stage的
            layers = Layer(dim=embed_dim,
                            depth=depths[i_layer],
                            num_heads=num_heads[i_layer],
                            input_resolution=[img_size[-2]//patch_size[-2], img_size[-1]//patch_size[-1]],
                            window_size=window_size,
                            mlp_ratio=self.mlp_ratio,
                            qkv_bias=qkv_bias,
                            drop=drop_rate,
                            attn_drop=attn_drop_rate,
                            drop_path=dpr[sum(depths[:i_layer]):sum(depths[:i_layer + 1])],
                            norm_layer=norm_layer,
                            use_checkpoint=use_checkpoint,
                            pre_norm=pre_norm
                            )
            self.layers.append(layers)

        self.pos_embed = nn.Parameter(torch.zeros(1, img_size[0]//patch_size[-2]*(img_size[1]//patch_size[-1]), embed_dim))
        nn.init.trunc_normal_(self.pos_embed, std=.02)

        self.apply(self._init_weights)

    def _init_weights(self, m):
        if isinstance(m, nn.Linear):
            nn.init.trunc_normal_(m.weight, std=.02)
            if isinstance(m, nn.Linear) and m.bias is not None:
                nn.init.constant_(m.bias, 0)
        elif isinstance(m, nn.LayerNorm):
            nn.init.constant_(m.bias, 0)
            nn.init.constant_(m.weight, 1.0)

    def forward(self, x):
        # x: [B, C, H, W]
        B, N, C = x.shape
        H, W = 16, 16  # 根据patch或窗口大小确定
        x = x.view(B, H, W, C)
        x = x.reshape(B, -1, C)
        T = 1
        # x, T, H, W = self.patch_embed(x)  # x:[B, H*W, C]
        x = x + self.pos_embed
        x = self.pos_drop(x)

        x = x.view(B, H, W, -1)

        for layer in self.layers:
            x= layer(x)

        return x

class LGUnet_all(nn.Module):
    """
    整个系统的主架构，集成了：
    编码器 Enc_net：局部特征提取；
    全局建模网络 LG_net：处理编码器提取的高层语义；
    解码器 Dec_net：输出最终预测图。
    """
    def __init__(self, img_size=[256, 256], patch_size=4, stride=[4,4], in_chans=4, out_chans=4, enc_depths=[2,2,6], enc_heads=[3,6,12], lg_depths=[2,2], lg_heads=[6,12], inchans_list=[4], outchans_list=[4], enc_dim=96,
                 embed_dim=768, window_size=8, Weather_T=16, drop_rate=0., attn_drop_rate=0., drop_path=0., use_checkpoint=False, channel_num=4, inp_length=1, use_mlp=False, pre_norm=True) -> None:
        super().__init__()

        # 局部编码器
        self.enc = Enc_net(img_size=img_size, patch_size=patch_size, stride=stride, inchans_list=inchans_list, embed_dim=enc_dim, lgnet_embed=embed_dim, depths=enc_depths,
                            num_heads=enc_heads, window_size=window_size, mlp_ratio=4, qkv_bias=True, drop_rate=drop_rate,
                            attn_drop_rate=attn_drop_rate, drop_path_rate=drop_path, use_checkpoint=use_checkpoint,
                            inp_length=inp_length, pre_norm=pre_norm)
        # 全局建模
        lg_patch_size = [stride[i] * 2** (len(enc_depths)-1) for i in range(len(stride))]
        self.net = LG_net(img_size=img_size, patch_size=lg_patch_size, \
                                        embed_dim=embed_dim, depths=lg_depths, num_heads=lg_heads, \
                                        window_size=window_size, drop_rate=drop_rate, attn_drop_rate=attn_drop_rate,
                                        drop_path_rate=drop_path, use_checkpoint=use_checkpoint, pre_norm=pre_norm)
        # 解码器
        self.dec = Dec_net(img_size=img_size, patch_size=patch_size, stride=stride, outchans_list=outchans_list,
                            channel_num=channel_num*2, embed_dim=enc_dim, lgnet_embed=embed_dim, depths=enc_depths,
                            num_heads=enc_heads, window_size=window_size, mlp_ratio=4, qkv_bias=True, drop_rate=drop_rate,
                            attn_drop_rate=attn_drop_rate, drop_path_rate=drop_path, use_checkpoint=use_checkpoint,
                            pre_norm=pre_norm)

        self.apply(self._init_weights)

    def _init_weights(self, m):
        if isinstance(m, nn.Linear):
            trunc_normal_(m.weight, std=.02)
            if isinstance(m, nn.Linear) and m.bias is not None:
                nn.init.constant_(m.bias, 0)
        elif isinstance(m, nn.LayerNorm):
            nn.init.constant_(m.bias, 0)
            nn.init.constant_(m.weight, 1.0)

    def forward(self, data):
        # data: [B, C, H, W]
        data = self.enc(data)     # 局部编码 -> [B, L, lgnet_embed]
        out = self.net(data)      # 全局建模 -> [B, H', W', C] 或视具体实现而定
        out = self.dec(out)       # 解码 -> [B, out_channels, H, W]
        return out

    def forward_get_latent(self, data):
        """
        返回解码结果与 latent（全局表征）。
        """
        data = self.enc(data)     # 编码
        latent = self.net(data)   # latent 表征
        out = self.dec(latent)    # 解码器输出
        return out, latent
    
    def encode(self, data):
        data = self.enc(data)
        latent = self.net(data)
        return latent

    def decode(self, latent):
        decoded = self.dec(latent)
        return decoded
