import os
import torch
import torch.nn as nn
from torch.utils.data import DataLoader
from torchinfo import summary
import torch.optim as optim
from dataset import WRFSlidingDataset
from vaelda.networks.transformer import LGUnet_all
from tqdm import tqdm
import matplotlib
matplotlib.use('Agg')  # ✅ 非交互式后端，适合服务器
import matplotlib.pyplot as plt
import numpy as np
from pytorch_msssim import ssim

import warnings
warnings.filterwarnings("ignore")

# === 组合损失 ===
def combined_loss(outputs, target, mask):
    outputs = outputs * mask
    target = target * mask
    mse_loss = nn.functional.mse_loss(outputs, target)
    ssim_loss = 1 - ssim(outputs, target, data_range=1.0, size_average=True)
    return mse_loss + 0.2 * ssim_loss

# === Accuracy 计算（可选，对回归类任务可定义阈值内为准确） ===
def calculate_accuracy(pred, target, threshold=0.1):
    correct = (torch.abs(pred - target) < threshold).float()
    return correct.mean().item()

def train_one_epoch(model, dataloader, optimizer, device):
    model.train()
    running_loss = 0.0
    for inputs, targets, mask in tqdm(dataloader, desc="Training"):
        inputs, targets, mask = inputs.to(device), targets.to(device), mask.to(device)

        optimizer.zero_grad()
        outputs = model(inputs)

        # print("outputs shape:", outputs.shape)
        # print("mask shape:", mask.shape)

        loss = combined_loss(outputs, targets, mask)
        loss.backward()
        optimizer.step()

        running_loss += loss.item() * inputs.size(0)

    return running_loss / len(dataloader.dataset)

@torch.no_grad()
def validate_one_epoch(model, dataloader, device):
    model.eval()
    total_loss, total_mae, total_rmse = 0, 0, 0
    total_acc, total_pixels = 0, 0

    for inputs, targets, mask in tqdm(dataloader, desc="Validation"):
        inputs, targets, mask = inputs.to(device), targets.to(device), mask.to(device)
        outputs = model(inputs)

        loss = nn.functional.mse_loss(outputs * mask, targets * mask)
        mae = nn.functional.l1_loss(outputs * mask, targets * mask)
        rmse = torch.sqrt(nn.functional.mse_loss(outputs * mask, targets * mask))

        acc = ((torch.abs(outputs - targets) < 0.5).float() * mask).sum()
        pixel_count = mask.sum()

        total_loss += loss.item() * inputs.size(0)
        total_mae += mae.item() * inputs.size(0)
        total_rmse += rmse.item() * inputs.size(0)
        total_acc += acc.item()
        total_pixels += pixel_count.item()

    avg_mse = total_loss / len(dataloader.dataset)
    avg_mae = total_mae / len(dataloader.dataset)
    avg_rmse = total_rmse / len(dataloader.dataset)
    avg_acc = total_acc / total_pixels

    return avg_mse, avg_mae, avg_rmse, avg_acc

def main():
    # === 配置参数 ===
    nc_dir = '/data/zsq_data/d01_d02_d03/d03_wrf'
    input_len = 1
    label_len = 1
    input_vars = ['U10', 'V10', 'T2', 'Q2']
    label_vars = ['U10', 'V10', 'T2', 'Q2']
    batch_size = 256
    epochs = 200
    lr = 1e-4

    # === 数据加载 ===
    dataset = WRFSlidingDataset(
        nc_dir=nc_dir,
        input_len=input_len,
        label_len=label_len,
        input_vars=input_vars,
        label_vars=label_vars
    )

    # === 时间顺序划分数据集（前10个季度训练，第11验证，第12测试）===
    quarter_ids = dataset.get_sample_quarters(base_year=2022)

    train_indices = [i for i, q in enumerate(quarter_ids) if 1 <= q <= 10]
    val_indices = [i for i, q in enumerate(quarter_ids) if q == 11]
    test_indices = [i for i, q in enumerate(quarter_ids) if q == 12]

    train_dataset = torch.utils.data.Subset(dataset, train_indices)
    val_dataset = torch.utils.data.Subset(dataset, val_indices)
    test_dataset = torch.utils.data.Subset(dataset, test_indices)

    print(f"[划分结果] 总样本数: {len(dataset)}")
    print(f"Train: {len(train_dataset)}  Val: {len(val_dataset)}  Test: {len(test_dataset)}")

    # === 构建 DataLoader ===
    train_loader = DataLoader(train_dataset, batch_size=batch_size, shuffle=True, num_workers=32, pin_memory=True)
    val_loader = DataLoader(val_dataset, batch_size=batch_size, shuffle=False, num_workers=32, pin_memory=True)
    test_loader = DataLoader(test_dataset, batch_size=batch_size, shuffle=False, num_workers=32, pin_memory=True)



    # === 模型构造 ===
    in_chans = len(input_vars) * input_len
    out_chans = len(label_vars) * label_len

    model = LGUnet_all(
        img_size=[256, 256],
        patch_size=4,
        stride=[4, 4],
        in_chans=in_chans,
        out_chans=out_chans,
        enc_depths=[2, 2, 6],
        enc_heads=[3, 6, 12],
        lg_depths=[2, 2],
        lg_heads=[6, 12],
        inchans_list=[in_chans],
        outchans_list=[out_chans],
        enc_dim=96,
        embed_dim=768,
        window_size=8,
        Weather_T=input_len,
        use_checkpoint=False,
        pre_norm=True
    )

    # === 多GPU设置 ===
    device_ids = [0, 1, 2, 3]  # 使用 GPU 0、1、2
    if torch.cuda.device_count() >= len(device_ids):
        print(f"✅ Using GPUs: {device_ids}")
        model = nn.DataParallel(model, device_ids=device_ids, output_device=device_ids[0])
        device = torch.device(f"cuda:{device_ids[0]}")
    else:
        device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

    model = model.to(device)

    # print("📌 模型结构如下：")
    # print(model)
    #
    # try:
    #     # 使用 torchinfo.summary 打印每层 shape（仅主GPU）
    #     summary(model, input_size=(batch_size, in_chans, 256, 256), device=device.type)
    # except Exception as e:
    #     print("⚠️ torchinfo.summary 显示失败，原因：", e)

    # === 损失函数 & 优化器 ===
    # criterion = nn.MSELoss()
    optimizer = torch.optim.Adam(model.parameters(), lr=lr)

    # 学习率调度器：监控 val_loss，如果 10 个 epoch 没有下降，就将 lr 乘以 0.5
    scheduler = optim.lr_scheduler.ReduceLROnPlateau(
        optimizer, mode='min',
        factor=0.5, patience=10
    )

    # === 初始化列表 ===
    train_losses, val_mse_list, val_mae_list, val_rmse_list, val_acc_list = [], [], [], [], []

    best_val_rmse = float('inf')

    # === 实时可视化设置 ===
    # plt.ion()
    fig, axs = plt.subplots(2, 2, figsize=(12, 8))
    axs = axs.ravel()
    axs[0].set_title('Train Loss')
    axs[1].set_title('Val Loss')
    axs[2].set_title('Val RMSE')
    axs[3].set_title('Val Accuracy')

    # === 训练主循环 ===
    for epoch in range(epochs):
        print(f"\nEpoch {epoch+1}/{epochs}")
        train_loss = train_one_epoch(model, train_loader, optimizer, device)
        val_mse, val_mae, val_rmse, val_acc = validate_one_epoch(model, val_loader, device)

        # 调度器根据验证集 MSE 更新学习率
        scheduler.step(val_mse)

        train_losses.append(train_loss)
        val_mse_list.append(val_mse)
        val_mae_list.append(val_mae)
        val_rmse_list.append(val_rmse)
        val_acc_list.append(val_acc)

        print(f"Train Loss: {train_loss:.4f} | Val MSE: {val_mse:.4f} | MAE: {val_mae:.4f} | RMSE: {val_rmse:.2f} | Acc: {val_acc:.4f}")

        # === 实时绘图 ===
        axs[0].plot(train_losses, label='Train Loss' if epoch == 0 else "")
        axs[1].plot(val_mae_list, label='Val MAE' if epoch == 0 else "")
        axs[2].plot(val_rmse_list, label='Val RMSE' if epoch == 0 else "")
        axs[3].plot(val_acc_list, label='Val Acc' if epoch == 0 else "")
        for ax in axs:
            ax.legend()
            ax.relim()
            ax.autoscale_view()
        plt.tight_layout()

        # 保存当前图像
        plt.savefig('metrics_forgoal.png')

        # === 保存最佳模型 ===
        if val_rmse < best_val_rmse:
            best_val_rmse = val_rmse
            torch.save(model.state_dict(), 'model_forgoal.pth')
            print(f"✅ Saved best model at epoch {epoch+1}")
    plt.close()

    print("\n🔍 Loading best model for final test...")
    model.load_state_dict(torch.load('best_model.pth'))
    model.eval()

    test_mse, test_mae, test_rmse, test_acc = validate_one_epoch(model, test_loader, device)
    print(f"\n🎯 Final Test Results: MSE: {test_mse:.4f} | MAE: {test_mae:.4f} | RMSE: {test_rmse:.2f} | Accuracy: {test_acc:.4f}")


if __name__ == "__main__":
    main()