from __future__ import annotations

import torch
import torch.nn as nn
import torch.nn.functional as F


def _pair_rotate(x: torch.Tensor) -> torch.Tensor:
    even = x[..., 0::2]
    odd = x[..., 1::2]
    return torch.stack((-odd, even), dim=-1).flatten(-2)


def _apply_rope_1d(x: torch.Tensor, pos: torch.Tensor) -> torch.Tensor:
    dim = x.shape[-1]
    if dim < 2:
        return x
    pairs = dim // 2
    freq = torch.arange(pairs, device=x.device, dtype=x.dtype)
    freq = 1.0 / (10000.0 ** (freq / max(pairs, 1)))
    theta = pos.to(dtype=x.dtype)[:, None] * freq[None, :]
    cos = torch.repeat_interleave(theta.cos(), 2, dim=-1)
    sin = torch.repeat_interleave(theta.sin(), 2, dim=-1)
    return x * cos[None, None, :, :] + _pair_rotate(x) * sin[None, None, :, :]


def apply_2d_rope(x: torch.Tensor, height: int, width: int) -> torch.Tensor:
    """Apply separable 2-D rotary embedding to `[B, heads, H*W, head_dim]`."""
    head_dim = x.shape[-1]
    if head_dim % 4 != 0:
        return x
    y, grid_x = torch.meshgrid(
        torch.arange(height, device=x.device),
        torch.arange(width, device=x.device),
        indexing="ij",
    )
    y = y.reshape(-1)
    grid_x = grid_x.reshape(-1)
    y_part, x_part = x.split(head_dim // 2, dim=-1)
    return torch.cat((_apply_rope_1d(y_part, y), _apply_rope_1d(x_part, grid_x)), dim=-1)


class CrossAttention2d(nn.Module):
    def __init__(self, channels: int, heads: int = 4, dropout: float = 0.0, backend: str = "auto") -> None:
        super().__init__()
        if channels % heads != 0:
            raise ValueError("channels must be divisible by heads")
        self.channels = channels
        self.heads = heads
        self.head_dim = channels // heads
        self.q = nn.Linear(channels, channels, bias=False)
        self.k = nn.Linear(channels, channels, bias=False)
        self.v = nn.Linear(channels, channels, bias=False)
        self.proj = nn.Linear(channels, channels)
        self.norm = nn.LayerNorm(channels)
        self.dropout = dropout
        self.backend = str(backend).lower()

    def _backend_order(self) -> tuple[str, ...]:
        if self.backend == "auto":
            return ("flash_attn", "xformers", "sdpa")
        return (self.backend,)

    def _attention(self, q: torch.Tensor, k: torch.Tensor, v: torch.Tensor) -> torch.Tensor:
        for backend in self._backend_order():
            if backend in {"flash", "flash_attn", "flash-attn"}:
                try:
                    from flash_attn import flash_attn_func

                    q_f = q.transpose(1, 2).contiguous()
                    k_f = k.transpose(1, 2).contiguous()
                    v_f = v.transpose(1, 2).contiguous()
                    out = flash_attn_func(
                        q_f,
                        k_f,
                        v_f,
                        dropout_p=self.dropout if self.training else 0.0,
                        causal=False,
                    )
                    return out.transpose(1, 2)
                except Exception:
                    continue
            if backend == "xformers":
                try:
                    import xformers.ops as xops

                    q_x = q.transpose(1, 2).contiguous()
                    k_x = k.transpose(1, 2).contiguous()
                    v_x = v.transpose(1, 2).contiguous()
                    out = xops.memory_efficient_attention(
                        q_x,
                        k_x,
                        v_x,
                        p=self.dropout if self.training else 0.0,
                    )
                    return out.transpose(1, 2)
                except Exception:
                    continue
            if backend in {"sdpa", "torch", "pytorch"}:
                return F.scaled_dot_product_attention(
                    q,
                    k,
                    v,
                    dropout_p=self.dropout if self.training else 0.0,
                    is_causal=False,
                )
        return F.scaled_dot_product_attention(
            q,
            k,
            v,
            dropout_p=self.dropout if self.training else 0.0,
            is_causal=False,
        )

    def forward(self, x: torch.Tensor, context: torch.Tensor) -> torch.Tensor:
        batch, channels, height, width = x.shape
        _, _, ctx_h, ctx_w = context.shape
        x_seq = x.permute(0, 2, 3, 1).reshape(batch, height * width, channels)
        ctx_seq = context.permute(0, 2, 3, 1).reshape(batch, ctx_h * ctx_w, channels)

        q = self.q(x_seq).view(batch, height * width, self.heads, self.head_dim).transpose(1, 2)
        k = self.k(ctx_seq).view(batch, ctx_h * ctx_w, self.heads, self.head_dim).transpose(1, 2)
        v = self.v(ctx_seq).view(batch, ctx_h * ctx_w, self.heads, self.head_dim).transpose(1, 2)
        q = apply_2d_rope(q, height, width)
        k = apply_2d_rope(k, ctx_h, ctx_w)

        out = self._attention(q, k, v)
        out = out.transpose(1, 2).reshape(batch, height * width, channels)
        out = self.norm(x_seq + self.proj(out))
        return out.reshape(batch, height, width, channels).permute(0, 3, 1, 2).contiguous()
