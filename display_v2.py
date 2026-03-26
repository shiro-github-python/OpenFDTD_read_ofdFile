import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button
from scipy.ndimage import zoom

def display_v2(v, display_res=300):
    """
    3次元numpy配列 v の断面を表示する関数。
    表示前に指定の画素数 (デフォルト300x300x300) に補完（リサイズ）します。
    v は v[x, y, z] の形式を想定しています。
    """
    # 1. 指定の画素数に合わせて補完 (Interpolation)
    current_shape = v.shape
    target_shape = (display_res, display_res, display_res)
    
    if current_shape != target_shape:
        zoom_factors = [t / c for t, c in zip(target_shape, current_shape)]
        print(f"Interpolating from {current_shape} to {target_shape}...")
        # 線形補完 (order=1) を使用
        v_zoomed = zoom(v, zoom_factors, order=1)
    else:
        v_zoomed = v

    # 2. MatplotlibによるUIの構築
    fig, ax = plt.subplots(figsize=(7, 8))
    plt.subplots_adjust(left=0.15, bottom=0.25)
    
    class State:
        def __init__(self):
            self.mode = 'XY'
            self.index = display_res // 2
            # ウィジェットの参照を保持（ガベージコレクション防止）
            self.slider = None
            self.btn_xy = None
            self.btn_yz = None
            self.btn_zx = None
            
    state = State()
    
    # 初回表示 (XY断面, z=index)
    # imshow(origin='lower') では data[y, x] の順になるため転置して表示
    im = ax.imshow(v_zoomed[:, :, state.index].T, origin='lower', cmap='viridis')
    ax.set_title(f"{state.mode} plane at index {state.index}")
    ax.set_xlabel('X')
    ax.set_ylabel('Y')

    # スライダの配置
    ax_slider = plt.axes([0.2, 0.12, 0.6, 0.03])
    state.slider = Slider(ax_slider, 'Index', 0, display_res - 1, valinit=state.index, valfmt='%d')

    def update_image():
        idx = int(state.index)
        if state.mode == 'XY':
            # 横軸: X, 縦軸: Y (固定: Z)
            data = v_zoomed[:, :, idx].T
            ax.set_xlabel('X')
            ax.set_ylabel('Y')
        elif state.mode == 'YZ':
            # Horizontal: Y, Vertical: Z. data[z, y] = v[x, y, z]
            data = v_zoomed[idx, :, :].T
            ax.set_xlabel('Y')
            ax.set_ylabel('Z')
        elif state.mode == 'ZX':
            # Horizontal: Z, Vertical: X. data[x, z] = v[x, y, z]
            data = v_zoomed[:, :, idx] # Wait, ZX plane usually means X and Z axes.
            # If we want Z horizontal and X vertical: data[x, z] = v[x, y, z]
            data = v_zoomed[:, idx, :] # Shape (X, Z). data[x, z] -> x vertical, z horizontal.
            ax.set_xlabel('Z')
            ax.set_ylabel('X')
        
        im.set_data(data)
        im.set_clim(vmin=np.min(data), vmax=np.max(data))
        ax.set_title(f"{state.mode} plane at index {idx}")
        fig.canvas.draw_idle()

    def on_slider_change(val):
        state.index = val
        update_image()

    state.slider.on_changed(on_slider_change)

    # Buttons
    ax_xy = plt.axes([0.2, 0.03, 0.15, 0.05])
    ax_yz = plt.axes([0.42, 0.03, 0.15, 0.05])
    ax_zx = plt.axes([0.65, 0.03, 0.15, 0.05])
    
    state.btn_xy = Button(ax_xy, 'XY')
    state.btn_yz = Button(ax_yz, 'YZ')
    state.btn_zx = Button(ax_zx, 'ZX')

    def set_xy(event):
        state.mode = 'XY'
        update_image()

    def set_yz(event):
        state.mode = 'YZ'
        update_image()

    def set_zx(event):
        state.mode = 'ZX'
        update_image()

    state.btn_xy.on_clicked(set_xy)
    state.btn_yz.on_clicked(set_yz)
    state.btn_zx.on_clicked(set_zx)

    plt.show()

# Example usage (can be removed or commented out)
if __name__ == "__main__":
    # Create dummy 1200x1200x1200 data if needed, but for testing let's use smaller
    # and let the function interpolate it.
    # Requirement example: 1200x1200x1200 -> 300x300x300
    print("Creating sample data...")
    x = np.linspace(0, 1, 120)
    y = np.linspace(0, 1, 120)
    z = np.linspace(0, 1, 120)
    X, Y, Z = np.meshgrid(x, y, z, indexing='ij')
    v_sample = np.sin(2 * np.pi * X) * np.cos(2 * np.pi * Y) * np.sin(2 * np.pi * Z)
    
    display_v(v_sample, display_res=30) # Smaller resolution for demonstration
