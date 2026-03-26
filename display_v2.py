import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button
from scipy.ndimage import zoom

def display_v2(v, display_res=300):
    """
    This function displays a cross-section of a 3D NumPy array `v`.
    Before displaying, it interpolates (resizes) to the specified number of pixels (default: 300x300x300).
    `v` is assumed to be in the format v[x, y, z].
    """
    # 1. Interpolate to match the specified number of pixels.(Interpolation)
    current_shape = v.shape
    target_shape = (display_res, display_res, display_res)
    
    if current_shape != target_shape:
        zoom_factors = [t / c for t, c in zip(target_shape, current_shape)]
        print(f"Interpolating from {current_shape} to {target_shape}...")
        # Use linear interpolation(order=1)
        v_zoomed = zoom(v, zoom_factors, order=1)
    else:
        v_zoomed = v

    # 2. Building a UI with Matplotlib
    fig, ax = plt.subplots(figsize=(8, 8))
    plt.subplots_adjust(left=0.15, bottom=0.25, right=0.85)
    
    class State:
        def __init__(self):
            self.mode = 'XY'
            self.index = display_res // 2
            # Retain widget references (prevent garbage collection)
            self.slider = None
            self.btn_xy = None
            self.btn_yz = None
            self.btn_zx = None
            
    state = State()
    
    # First display (XY cross section, z=index)
    # Since imshow(origin='lower') displays the data in the order data[y, x], it needs to be transposed before display.
    im = ax.imshow(v_zoomed[:, :, state.index].T, origin='lower', cmap='viridis')
    ax.set_title(f"{state.mode} plane at index {state.index}")
    ax.set_xlabel('X')
    ax.set_ylabel('Y')

    # colorbar placement
    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label('Value')

    # Slider placement
    ax_slider = plt.axes([0.2, 0.12, 0.6, 0.03])
    state.slider = Slider(ax_slider, 'Index', 0, display_res - 1, valinit=state.index, valfmt='%d')

    def update_image():
        idx = int(state.index)
        if state.mode == 'XY':
            # Horizontal: X, Vertical: Y. data[y, x] = v[x, y, z]
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
            data = v_zoomed[:, idx, :]
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

