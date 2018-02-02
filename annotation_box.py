# from https://matplotlib.org/examples/pylab_examples/demo_annotation_box.html
import matplotlib.pyplot as plt
import numpy as np

from matplotlib.offsetbox import OffsetImage, AnnotationBbox

fig, ax = plt.subplots()

xy = (0.5, 0.5)
arr = np.arange(100).reshape((10, 10))
im = OffsetImage(arr, zoom=5)
im.image.axes = ax

ab = AnnotationBbox(im, xy,
                    xybox=(-50., 70.),
                    xycoords='data',
                    boxcoords="offset points",
                    pad=0.3,
                    arrowprops=dict(arrowstyle="->"))

ax.add_artist(ab)

ax.set_xlim(0, 1)
ax.set_ylim(0, 1)

plt.show()
