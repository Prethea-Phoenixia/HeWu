import matplotlib as mpt
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib.ticker import MultipleLocator
from math import floor, log10, sqrt, pi
import numpy as np
import numpy.ma as ma
from HeWu.modelWE1984 import therm, phi


def secant(f, x_0, x_1, x_min=None, x_max=None, tol=1e-6, it=1000):
    """secant method that solves f(x) = 0 subjected to x in [x_min,x_max]"""
    if x_min is not None:
        if x_0 < x_min:
            x_0 = x_min
        if x_1 < x_min:
            x_1 = x_min
    if x_max is not None:
        if x_0 > x_max:
            x_0 = x_max
        if x_1 > x_max:
            x_1 = x_max

    fx_0 = f(x_0)
    fx_1 = f(x_1)
    for _ in range(it):
        x_2 = x_1 - fx_1 * (x_1 - x_0) / (fx_1 - fx_0)
        if x_min is not None and x_2 < x_min:
            x_2 = x_min
        if x_max is not None and x_2 > x_max:
            x_2 = x_max
        x_0, x_1, fx_0, fx_1 = x_1, x_2, fx_1, f(x_2)
        if abs(fx_1) < tol:
            return x_1, fx_1

    raise ValueError("Maximum iteration exceeded at ({},{})".format(x_1, fx_1))


def qx(x):
    return therm(Y, 0, x, vis)


def qy(y):
    return therm(Y, y, 0, vis)


clevels = (
    "gold",
    "gold",
    "orange",
    "orange",
    "red",
    "red",
    "darkred",
)


vis = 15000  # 10km -80km
Y = 1
Y3 = Y ** (1 / 3)
logY = log10(Y)

burns = (
    0.26 * logY + 1.1,
    0.41 * logY + 2,
    0.60 * logY + 2.85,
    0.75 * logY + 4,
    0.88 * logY + 5,
    1.25 * logY + 6,
    1.45 * logY + 7,
)


ymax, _ = secant(lambda y: qy(y) - (0.26 * logY + 1), 1, 2, x_min=0)
xmax = ymax

delta = 2 * Y3
rows = floor(ymax // delta) + 1
columns = floor(xmax // delta) + 1
x = delta * np.array(range(columns)) / 1e3
y = delta * np.array(range(rows)) / 1e3
F = np.zeros((rows, columns))

prob1 = []
mu1 = burns[2]
sigma1L = (mu1 - burns[1]) / 1.17741
sigma1H = (-mu1 + burns[3]) / 1.17741

prob2 = []
mu2 = burns[4]
sigma2L = (mu2 - burns[3]) / 1.17741
sigma2H = (-mu2 + burns[5]) / 1.17741

prob3 = []
mu3 = burns[6]
sigma3 = (mu3 - burns[5]) / 1.17741

for i in range(columns):
    groundRange = i * delta
    for j in range(rows):
        height = j * delta
        if i == 0 and j == 0:
            groundRange += 0.1
            height += 0.1

        flu = therm(Y, height, groundRange, vis)
        F[j][i] = flu

        if i == 0:
            if flu < mu1:
                prob1.append(phi((flu - mu1) / sigma1L) * sqrt(2 * pi))
            else:
                prob1.append(phi((flu - mu1) / sigma1H) * sqrt(2 * pi))
            if flu < mu2:
                prob2.append(phi((flu - mu2) / sigma2L) * sqrt(2 * pi))
            else:
                prob2.append(phi((flu - mu2) / sigma2H) * sqrt(2 * pi))
            if flu < mu3:
                prob3.append(phi((flu - mu3) / sigma3) * sqrt(2 * pi))
            else:
                prob3.append(1)


fig, ax = plt.subplots(1, 1, figsize=(11.7, 8.3))

rect = patches.Rectangle(
    (0, 0), 2.2 * Y3, 1.5 * Y3, linewidth=1.5, edgecolor="r", facecolor="none"
)
ax.add_patch(rect)
ax.text(
    2.2 * Y3,
    1.5 * Y3,
    "Model Limits",
    color="red",
    va="top",
    ha="right",
    alpha=0.75,
)

i = 0
for b in burns:
    gr, _ = secant(lambda gr: qx(gr) - b, 1, 2)
    ax.plot(gr / 1000, 0, marker=7, color=clevels[i], markersize=8, zorder=4)
    i += 1

csf = ax.contour(
    x,
    y,
    F,
    colors=clevels,
    levels=burns,
    linestyles=(
        "solid",
        "dotted",
    ),
)


csh = ax.contourf(x, y, F, colors=clevels, levels=burns, extend="max", alpha=0.33)

i = 0
burnh = []
for b in burns:
    h, _ = secant(lambda y: qy(y) - b, 1, 2, x_min=0)
    ax.text(
        1 * Y3 / 1e3,
        (h - 10 * Y3) / 1e3,
        "{:.1f}".format(burns[i]),
        va="top",
        ha="left",
        alpha=0.75,
    )
    burnh.append(h)
    i += 1


ax.set_aspect("equal")
ax.grid(which="major", color="#DDDDDD", linewidth=1)
ax.grid(which="minor", color="#EEEEEE", linestyle=":", linewidth=1)
ax.minorticks_on()


nmf, lblf = csf.legend_elements()
lblf = (
    "<1°",
    "<1°-1°",
    "1°",
    "1°-2°",
    "2°",
    "2°-3°",
    "3°",
)

legend2 = ax.legend(
    nmf[::], lblf, title="Burn to Exposed Personnel", ncol=1, loc="lower right"
)
ax.set_title("Thermal Fluence (Cal/cm$^2$)")

ax.set_xlabel("Ground Range (km)")
fig.suptitle(
    "{:.0f}kT Airburst -{:.0f}km Visibility".format(Y, vis / 1000), fontsize="x-large"
)


_, top = ax.get_ylim()

axins = inset_axes(
    ax,
    width="100%",
    height="100%",
    loc="lower right",
    borderpad=0,
    bbox_to_anchor=(-1.25 * top / 10, 0, top / 10, top),
    bbox_transform=ax.transData,
)
probs = (prob1, prob2, prob3)
for i in range(3):
    axins.plot(
        (0, 0.5),
        (burnh[2 * i + 1], burnh[2 * i + 1]),
        c=clevels[2 * i + 1],
        linestyle="dotted",
        linewidth=2,
    )
    axins.plot(probs[i], y, c=clevels[2 * i + 2])
    axins.fill_betweenx(y, probs[i], 0, color=clevels[2 * i + 2], alpha=0.2)

axins.set_xlim(left=0, right=1)
axins.set_ylim(bottom=0, top=top)
axins.set_xlabel("Burn Probability")
axins.minorticks_on()
axins.grid(which="major", color="#DDDDDD", linewidth=1)
axins.grid(which="minor", color="#EEEEEE", linestyle=":", linewidth=1)
axins.xaxis.set_minor_locator(MultipleLocator(0.25))
axins.invert_xaxis()
ax.set_yticklabels([])
axins.set_ylabel("Burst Height (km)")
plt.show()
