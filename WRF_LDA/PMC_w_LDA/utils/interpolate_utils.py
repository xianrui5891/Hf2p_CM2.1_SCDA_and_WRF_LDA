import torch

def bilinear_interpolate_torch(corners, values, pt, maxiter=20, tol=1e-8):
    """
    简化版：如果 grid 是规则轴向的，直接用 fractional s,t 更简单。
    但这里保留接口：如果你坚持用逆映射求 s,t，这里用 Newton 的 torch 版本（逐点）。
    corners: tensor (4,2)
    values: tensor (4,)
    pt: tensor (2,)
    返回: scalar tensor (插值值)
    """
    # Newton 求解 (s,t) 使得 map_st(s,t) = pt
    # 初值
    s = torch.tensor(0.5, device=pt.device, dtype=pt.dtype, requires_grad=False)
    t = torch.tensor(0.5, device=pt.device, dtype=pt.dtype, requires_grad=False)

    p0, p1, p2, p3 = corners[0], corners[1], corners[2], corners[3]

    def map_st(sv, tv):
        return (1 - sv) * (1 - tv) * p0 + sv * (1 - tv) * p1 + sv * tv * p2 + (1 - sv) * tv * p3

    def jacobian(sv, tv):
        # 返回 2x2 tensor
        dx_ds = (-(1 - tv)) * p0 + (1 - tv) * p1 + tv * p2 + (-tv) * p3
        dx_dt = (-(1 - sv)) * p0 + (-sv) * p1 + sv * p2 + (1 - sv) * p3
        J = torch.stack([dx_ds, dx_dt], dim=1)  # (2,2)
        return J

    for i in range(maxiter):
        F = map_st(s, t) - pt  # (2,)
        if torch.norm(F) < tol:
            break
        J = jacobian(s, t)
        # solve J delta = -F
        try:
            delta = torch.linalg.solve(J, -F)  # (2,)
        except RuntimeError:
            break
        s = s + delta[0]; t = t + delta[1]
    # 权重
    w0 = (1 - s) * (1 - t)
    w1 = s * (1 - t)
    w2 = s * t
    w3 = (1 - s) * t
    val = w0 * values[0] + w1 * values[1] + w2 * values[2] + w3 * values[3]
    return val

def mean_value_coordinates_interp_torch(corners, values, pt, eps=1e-15):
    """
    PyTorch 版 mean-value coordinates for polygon (2D)
    corners: (n,2), values: (n,), pt: (2,)
    返回 scalar tensor
    """
    p = corners  # (n,2)
    v = values   # (n,)
    x = pt
    r = p - x  # (n,2)
    d = torch.linalg.norm(r, dim=1)  # (n,)
    # 如果靠近顶点
    idx_zero = torch.where(d < eps)[0]
    if idx_zero.numel() > 0:
        return v[int(idx_zero[0])]
    n = p.shape[0]
    # 计算 theta_i between r_i and r_{i+1}
    thetas = []
    for i in range(n):
        u = r[i]
        wv = r[(i + 1) % n]
        cuw = torch.dot(u, wv)
        suw = u[0]*wv[1] - u[1]*wv[0]  # cross scalar in 2D
        ang = torch.atan2(torch.abs(suw), cuw)
        thetas.append(ang)
    thetas = torch.stack(thetas)  # (n,)
    w = torch.zeros(n, dtype=p.dtype, device=p.device)
    for i in range(n):
        im = (i - 1) % n
        denom = d[i]
        if denom < eps:
            w[i] = 1e18
        else:
            w[i] = (torch.tan(thetas[im] / 2.0) + torch.tan(thetas[i] / 2.0)) / denom
    S = torch.sum(w)
    if S == 0:
        return torch.tensor(float('nan'), device=p.device, dtype=p.dtype)
    lambda_w = w / S
    val = torch.dot(lambda_w, v)
    return val

def interpolate(val_field: torch.Tensor, axis: list, coors: list, type: str = 'bilinear'):
    """
    torch 版 interpolate，接口尽量与 numpy 版相同。
    val_field: torch tensor shape (nx, ny) or [i,j]
    axis: [axis_x (1D numpy-like), axis_y (1D numpy-like)]  这里允许为 numpy 或 torch
    coors: iterable of (x,y) pairs (can be lists, tuples or torch tensors)
    返回: list of scalar tensors -> 最终你可以 torch.stack(...) 变成 1D tensor
    注意：不做 detach，不做 CPU 转换，全部保持在 val_field.device 上。
    """
    device = val_field.device
    dtype = val_field.dtype

    # 确保 axis 是 torch tensor
    ax0 = torch.as_tensor(axis[0], device=device, dtype=dtype)
    ax1 = torch.as_tensor(axis[1], device=device, dtype=dtype)

    nx = ax0.numel(); ny = ax1.numel()

    values = []
    for (x_raw, y_raw) in coors:
        x = float(x_raw) if not torch.is_tensor(x_raw) else x_raw.item() if x_raw.numel()==1 else x_raw
        y = float(y_raw) if not torch.is_tensor(y_raw) else y_raw.item() if y_raw.numel()==1 else y_raw

        # 只处理落在内部的点（和你原代码一致）
        if (ax0[0] < x < ax0[-1]) and (ax1[0] < y < ax1[-1]):
            # 找到 idx, idy: searchsorted 等价
            # torch.searchsorted 需要 1D ascending tensor
            idx = int(torch.searchsorted(ax0, torch.tensor(x, device=device, dtype=dtype), right=True).item())
            idy = int(torch.searchsorted(ax1, torch.tensor(y, device=device, dtype=dtype), right=True).item())

            print(f"the interpolate indexes are: idx: {idx}, idy: {idy}.")

            # 边界保护（与 numpy 代码一致）
            if idx == 0 or idy == 0 or idx >= nx or idy >= ny:
                continue

            # corner coordinates as in your numpy code
            corners = torch.stack([
                torch.stack([ax0[idx - 1], ax1[idy - 1]]),
                torch.stack([ax0[idx],     ax1[idy - 1]]),
                torch.stack([ax0[idx - 1], ax1[idy]]),
                torch.stack([ax0[idx],     ax1[idy]])
            ], dim=0)  # (4,2)

            corner_values = torch.stack([
                val_field[idx - 1, idy - 1],
                val_field[idx    , idy - 1],
                val_field[idx - 1, idy    ],
                val_field[idx    , idy    ]
            ])  # (4,)

            pt = torch.tensor([x, y], device=device, dtype=dtype)

            if type == 'bilinear':
                # 更简单且常用：直接用 fractional s,t：
                x0 = ax0[idx - 1].item(); x1 = ax0[idx].item()
                y0 = ax1[idy - 1].item(); y1 = ax1[idy].item()
                # 注意：把 s,t 定义为 floats，但插值权重在 torch 下算（保持 val_field 的 device）
                # 为了兼容任意轴间距，我们用 fractions：
                denom_x = (x1 - x0) if (x1 - x0) != 0 else 1.0
                denom_y = (y1 - y0) if (y1 - y0) != 0 else 1.0
                s = (x - x0) / denom_x
                t = (y - y0) / denom_y
                # 转为 torch scalar 保持 device/dtype
                s_t = torch.as_tensor([s, t], device=device, dtype=dtype)
                s = s_t[0]; t = s_t[1]
                w0 = (1 - s) * (1 - t)
                w1 = s * (1 - t)
                w2 = (1 - s) * t
                w3 = s * t
                val = w0 * corner_values[0] + w1 * corner_values[1] + w2 * corner_values[2] + w3 * corner_values[3]
            elif type == 'mean_value':
                val = mean_value_coordinates_interp_torch(corners, corner_values, pt)
            else:
                raise RuntimeError("this type is not supported!")
            values.append(val)
    return values  # list of scalar torch tensors