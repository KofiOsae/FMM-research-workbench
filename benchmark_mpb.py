"""Compare the 2D band solver with the published MPB square-rod example."""

import numpy as np
from scipy.linalg import eigh

from bands import BandModel, _convolution, dielectric_grid


def mpb_te_first_band(fourier_order: int, grid_size: int = 512) -> float:
    model = BandModel(inclusion_n=np.sqrt(12), fill=.2, fourier_order=fourier_order,
                      grid_size=grid_size)
    eps = dielectric_grid(model)
    orders = np.arange(-fourier_order, fourier_order+1)
    mx, my = np.meshgrid(orders, orders, indexing="ij")
    mx, my = mx.ravel(), my.ravel()
    inverse_conv = _convolution(1/eps, mx, my)
    gx, gy = mx+.3, my+.3
    matrix = gx[:,None]*inverse_conv*gx[None,:] + gy[:,None]*inverse_conv*gy[None,:]
    return float(np.sqrt(max(eigh(matrix, eigvals_only=True, subset_by_index=(0,0))[0], 0)))


if __name__ == "__main__":
    for order in (3, 5, 7):
        print(order, mpb_te_first_band(order), "MPB reference: 0.372604")
