"""
debug_utils.py

提供一个简单的绘图类 DebugPlotter，用于在训练/迭代过程中记录并绘制三个 loss 随迭代次数的变化。

接口：
- plot(it, loss_J, loss_background, loss_observation)
    每次调用会把传入的数值追加到内部缓存，并在当前图上绘制五个子图：
      1) 三条 loss 的线性刻度组合图（带数据点），x 轴从 0 开始（避免出现 offset 科学计数显示为 0*1e7）
      2) 三条 loss 的对数刻度组合图（带数据点）
      3) loss_J 的单独变化图（带数据点）
      4) loss_background 的单独变化图（带数据点）
      5) loss_observation 的单独变化图（带数据点）

    绘图同时会把当前会话的数据保存为同一目录下的 npz 文件（文件名与图片相对应），
    当前会话未 close 之前会覆盖该 npz 与 png（保持最新），close() 后索引自增，下一次会产生新文件。

- close()
    关闭当前图像并清空内部缓存，下一次 plot 会新建文件并重新开始记录。

实现尽量简洁，不做过多鲁棒性检查。
"""

from typing import Optional, Tuple
import numpy as np
import matplotlib.pyplot as plt
import os



class DebugPlotter:
    """简单的绘图器。

    使用示例：
        dp = DebugPlotter('training_plot')
        dp.plot(it, loss_J, loss_bg, loss_obs)
        dp.close()

    注意：传入的 loss 参数应为标量（例如 float 或 tensor.item() 的返回值）。
    """

    def __init__(self, filename_prefix: str = 'training_plot', dpi: int = 150, figsize=(8, 10)):
        self.filename_prefix = filename_prefix
        self.dpi = dpi
        self.figsize = figsize

        # 输出目录（在当前工作目录下新建子文件夹，用于存储图片和数据）
        self.out_dir = f"{self.filename_prefix}_outputs"
        os.makedirs(self.out_dir, exist_ok=True)


        # 数据缓存
        self.iters = []
        self.loss_J = []
        self.loss_bg = []
        self.loss_obs = []

        # matplotlib 对象
        self._fig: Optional[plt.Figure] = None
        self._axs = None

        # 文件索引：close 后自增，保证新会话使用新文件
        self._file_index = 0

    def _open_figure_if_needed(self):
        if self._fig is None:
            # 5 行 1 列子图：前两个为组合图（线性/对数），后三个分别为单独 loss
            self._fig, axs = plt.subplots(5, 1, figsize=self.figsize)
            self._axs = axs
            self._fig.tight_layout()

    def plot(self, it, loss_J, loss_background, loss_observation):
        """追加数据并绘制（保存到文件）。

        返回 (image_filename, data_filename)
        """
        # 强制转换为 Python float（如果是 numpy / torch scalar）以避免绘图问题
        x = float(it)
        yJ = float(loss_J)
        ybg = float(loss_background)
        yobs = float(loss_observation)

        # 追加
        self.iters.append(x)
        self.loss_J.append(yJ)
        self.loss_bg.append(ybg)
        self.loss_obs.append(yobs)

    def close(self):
        # 打开（如尚未打开）
        self._open_figure_if_needed()

        # 清空并重画所有子图
        ax_comb_lin = self._axs[0]
        ax_comb_log = self._axs[1]
        ax_J = self._axs[2]
        ax_bg = self._axs[3]
        ax_obs = self._axs[4]

        for ax in self._axs:
            ax.cla()

        # 组合线性图（从 x=0 开始，禁止 offset 科学计数显示）
        ax_comb_lin.plot(self.iters, self.loss_J, marker='o', markersize=4, linewidth=1, label='loss_J')
        ax_comb_lin.plot(self.iters, self.loss_bg, marker='o', markersize=4, linewidth=1, label='loss_background')
        ax_comb_lin.plot(self.iters, self.loss_obs, marker='o', markersize=4, linewidth=1, label='loss_observation')
        ax_comb_lin.set_xlabel('iteration')
        ax_comb_lin.set_ylabel('loss')
        ax_comb_lin.legend()
        ax_comb_lin.grid(True, linestyle=':', linewidth=0.5)
        # 确保 x 从 0 开始，若只有一个点则给出合理范围
        x_max = max(self.iters) if len(self.iters) > 0 else 1
        ax_comb_lin.set_xlim(left=0, right=max(x_max, 1))
        ax_comb_lin.ticklabel_format(useOffset=False, style='plain')

        # 组合对数图
        ax_comb_log.plot(self.iters, self.loss_J, marker='o', markersize=4, linewidth=1, label='loss_J')
        ax_comb_log.plot(self.iters, self.loss_bg, marker='o', markersize=4, linewidth=1, label='loss_background')
        ax_comb_log.plot(self.iters, self.loss_obs, marker='o', markersize=4, linewidth=1, label='loss_observation')
        ax_comb_log.set_xlabel('iteration')
        ax_comb_log.set_ylabel('loss (log scale)')
        ax_comb_log.set_yscale('log')
        ax_comb_log.legend()
        ax_comb_log.grid(True, which='both', linestyle=':', linewidth=0.5)
        ax_comb_log.set_xlim(left=0, right=max(x_max, 1))

        # 单独三个子图（带点）
        ax_J.plot(self.iters, self.loss_J, marker='o', markersize=4, linewidth=1)
        ax_J.set_ylabel('loss_J')
        ax_J.set_xlim(left=0, right=max(x_max, 1))
        ax_J.grid(True, linestyle=':', linewidth=0.5)

        ax_bg.plot(self.iters, self.loss_bg, marker='o', markersize=4, linewidth=1)
        ax_bg.set_ylabel('loss_background')
        ax_bg.set_xlim(left=0, right=max(x_max, 1))
        ax_bg.grid(True, linestyle=':', linewidth=0.5)

        ax_obs.plot(self.iters, self.loss_obs, marker='o', markersize=4, linewidth=1)
        ax_obs.set_ylabel('loss_observation')
        ax_obs.set_xlabel('iteration')
        ax_obs.set_xlim(left=0, right=max(x_max, 1))
        ax_obs.grid(True, linestyle=':', linewidth=0.5)

        # 保存到文件（图片和数据）
        img_filename = os.path.join(self.out_dir, f"{self.filename_prefix}_{self._file_index}.png")
        self._fig.savefig(img_filename, dpi=self.dpi, bbox_inches='tight')

        # 将数据保存为 npz，覆盖当前会话文件；close() 后索引自增，下一会话新文件
        """关闭当前图像并清空内部缓存，下一次 plot 会新建文件并重新开始记录。"""
        if self._fig is not None:
            plt.close(self._fig)
            self._fig = None
            self._axs = None
        # 清空缓存，保证下一次为全新会话
        self.iters = []
        self.loss_J = []
        self.loss_bg = []
        self.loss_obs = []
        self._file_index += 1


# 简单测试（仅在直接运行此文件时执行，不会在导入时运行）
if __name__ == '__main__':
    dp = DebugPlotter()
    for i in range(1, 11):
        dp.plot(i, i * 0.5 + 1, 1.0 / i + 0.1, (i % 3) + 0.2)
    dp.close()
