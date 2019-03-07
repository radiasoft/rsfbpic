# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
import matplotlib.colors as colors

class RsMidpointNormalize(colors.Normalize):
    """
    Normalise the colorbar for 2D heatmap/contour plots.

    Colors will be used symmetrically around a prescribed midpoint value.
    Example of usage:
        im=ax1.imshow(array,norm=MidpointNormalize(midpoint=0.,vmin=-10,vmax=10))
    Args:
        colors.Normalize:    ...not sure (DLB)
    Returns:
        np.ma.masked_array:  ...not sure (DLB)
    """
    def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
        self.midpoint = midpoint
        colors.Normalize.__init__(self, vmin, vmax, clip)

    def __call__(self, value, clip=None):
        # ignore masked values and various edge cases for now...
        x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
        return np.ma.masked_array(np.interp(value, x, y), np.isnan(value))
