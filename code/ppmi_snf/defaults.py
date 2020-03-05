# -*- coding: utf-8 -*-
"""
Some colors that are commonly used for plotting figures throughout the project
"""

from matplotlib.colors import ListedColormap
import seaborn as sns

# colors and colormaps
gray = [
    (0.65098039, 0.65882352, 0.67058823, 1.0)
]
edgegray = [
    (0.24313725, 0.24313725, 0.24313725, 1.0)
]
similarity_cmap = [
    (0.94901961, 0.94901961, 0.94901961, 1.0),
    (0.50196078, 0.16470588, 1.0, 0.73725490),
]
three_cluster_cmap = [
    (0.14901960, 0.46274509, 0.72156862, 1.0),
    (0.15686274, 0.56470588, 0.28627450, 1.0),
    (0.86666666, 0.30196078, 0.01568627, 1.0)
]
four_cluster_cmap = three_cluster_cmap + [
    (1.0, 0.8, 0.0, 1.0)
]
similarity_cmap = ListedColormap(sns.blend_palette(similarity_cmap, 1000))

del ListedColormap, sns
