import matplotlib.pyplot as plt
import matplotlib as mpt
from math import floor, log10, sqrt, pi
import matplotlib.patches as patches
import numpy as np
import numpy.ma as ma
from HeWu.modelWE1984 import crater
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib.ticker import MultipleLocator


Y = 1000
Y3 = Y ** (1 / 3.1)


ymax = 3 * Y3
ymin = -75 * Y3

delta = 0.5 * Y3
ylevels = floor((ymax - ymin) // delta) + 1


fig, axd = plt.subplot_mosaic(
    [["upper left", "right"], ["lower left", "right"]],
    layout="constrained",
)
# ax2 = ax1.twiny()
ax1 = axd["right"]
ax2 = axd["upper left"]
ax3 = axd["lower left"]


# fig, axs = plt.subplots(1, 3, figsize=(11.7, 8.3), sharey=True)
# ax1, ax2, ax3 = axs
# ax1, ax2 = axs
H = []
for i in range(ylevels):
    height = ymax - i * delta
    H.append(height)

colors = ("brown", "dimgrey", "black")
linestyles = ("solid", "dotted")
labels = ("Dry Dirt", "Wet Dirt", "Dry Soft Rock", "Wet Soft Rock", "Hard Rock")

CVS = []
CRS = []
CDS = []
for ground in range(1, 6):
    CV = []
    CR = []
    CD = []
    for height in H:
        try:
            v, r, d, _, _, _ = crater(Y, height, ground, 1, ground, 1, ground, GR=None)
        except TypeError:
            v, r, d = 0, 0, 0
        CV.append(v / 1e3)
        CR.append(r)
        CD.append(d)

    CVS.append(CV)
    CRS.append(CR)
    CDS.append(CD)

    ax1.plot(
        CV,
        H,
        color=colors[(ground - 1) // 2],
        linestyle=linestyles[(ground - 1) % 2],
        label=labels[ground - 1],
    )
    ax2.plot(
        CR,
        H,
        color=colors[(ground - 1) // 2],
        linestyle=linestyles[(ground - 1) % 2],
        label=labels[ground - 1],
    )

    ax3.plot(
        CD,
        H,
        color=colors[(ground - 1) // 2],
        linestyle=linestyles[(ground - 1) % 2],
        label=labels[ground - 1],
    )


ax1.fill_betweenx(H, CVS[2], CVS[3], color=colors[1], alpha=0.1)
ax1.fill_betweenx(H, CVS[0], CVS[1], color=colors[0], alpha=0.1)
ax2.fill_betweenx(H, CRS[2], CRS[3], color=colors[1], alpha=0.1)
ax2.fill_betweenx(H, CRS[0], CRS[1], color=colors[0], alpha=0.1)
ax3.fill_betweenx(H, CDS[2], CDS[3], color=colors[1], alpha=0.1)
ax3.fill_betweenx(H, CDS[0], CDS[1], color=colors[0], alpha=0.1)
legend = ax1.legend(title="Geological Condition", ncol=1, loc="lower right")

ax1.set_ylim(ymin, ymax)
ax1.set_xlim(0)
_, xmax = ax1.get_xlim()
rect = patches.Rectangle(
    (0, -40 * Y3),
    xmax * 0.999,
    43 * Y3,
    linewidth=1.5,
    edgecolor="r",
    facecolor="none",
)
ax1.add_patch(rect)
ax1.minorticks_on()
ax1.grid(which="major", color="#DDDDDD", linewidth=1)
ax1.grid(which="minor", color="#EEEEEE", linestyle=":", linewidth=1)
ax1.set_xlabel("Apparent Crater Volume ($\\times 10^3$ m$^3$)")
ax1.set_ylabel("Burst Height (m)")

ax1.yaxis.tick_right()
ax1.yaxis.set_label_position("right")

ax1.text(
    xmax,
    ymax * 0.9,
    "Model Limits",
    va="top",
    ha="right",
    color="red",
    alpha=0.75,
)

ax2.set_xlabel("Apparent Crater Radius (m)")
ax2.xaxis.tick_top()
ax2.xaxis.set_label_position("top")
ax2.set_ylim(ymin, ymax)
ax2.tick_params(
    axis="x",  # changes apply to the x-axis
    which="both",  # both major and minor ticks are affected
    bottom=False,  # ticks along the bottom edge are off
    top=True,  # ticks along the top edge are off
    labelbottom=False,
)  # labels along the bottom edge are off
ax2.set_ylabel("Burst Height (m)")
ax2.grid(which="major", color="#DDDDDD", linewidth=1)
ax2.grid(which="minor", color="#EEEEEE", linestyle=":", linewidth=1)
ax2.set_xlim(0)
xmin, xmax = ax2.get_xlim()
rect = patches.Rectangle(
    (0, -40 * Y3),
    xmax,
    43 * Y3,
    linewidth=1.5,
    edgecolor="r",
    facecolor="none",
)
ax2.add_patch(rect)

ax3.set_xlim(0, xmax)
ax3.set_ylim(ymin, ymax)
ax3.set_xlabel("Apparent Crater Depth (m)")
ax3.grid(which="major", color="#DDDDDD", linewidth=1)
ax3.grid(which="minor", color="#EEEEEE", linestyle=":", linewidth=1)
ax3.set_ylabel("Burst Height (m)")
rect = patches.Rectangle(
    (0, -40 * Y3),
    xmax,
    43 * Y3,
    linewidth=1.5,
    edgecolor="r",
    facecolor="none",
)
ax3.add_patch(rect)


fig.suptitle("{:.0f}kT Underground Burst".format(Y), fontsize="x-large")
plt.show()
