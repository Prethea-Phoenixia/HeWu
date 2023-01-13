import numpy as np
import numpy.ma as ma
from math import floor
import matplotlib as mpt
from matplotlib import ticker
import matplotlib.pyplot as plt
from HeWu.modelBrode1987 import airburst
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib.patches as patches
from math import sqrt, log, log10


Y = 1
Y3 = Y ** (1 / 3)
xmax = 3400 * Y3
ymax = 3400 * Y3
xmed = 1200 * Y3
ymed = 1200 * Y3
delta = 10 * Y3

intervals = (100, 200, 250, 500, 1000, 2000, 2500, 5000, 10000)

for imed in intervals:
    if xmed / imed <= 10:
        break
for imax in intervals:
    if xmax / imax <= 10:
        break

lowXs = floor(xmed / delta)
lowYs = floor(ymed / delta)

x = np.append(
    np.arange(delta / 1e3, xmed / 1e3, delta / 1e3),
    np.arange(xmed / 1e3, (xmax + delta) / 1e3, 5 * delta / 1e3),
)  # in kilometer

y = np.append(
    np.arange(delta / 1e3, ymed / 1e3, delta / 1e3),
    np.arange(ymed / 1e3, (ymax + delta) / 1e3, 5 * delta / 1e3),
)  # in kilometer

rows = y.size
columns = x.size
P = np.zeros((rows, columns))
Q = np.zeros((rows, columns))  # dynamic pressure
T = np.zeros((rows, columns))  # overpressure arrival time
DPQ = np.zeros((rows, columns))  # dynamic pressure positive phase duration
DPP = np.zeros((rows, columns))  # overpressure positive phase duration
IP = np.zeros((rows, columns))  # overpressure total impulse
IQ = np.zeros((rows, columns))  # dynamic pressure total impulse

i = 0
for gr in x:
    groundRange = gr * 1e3
    j = 0
    for h in y:
        height = h * 1e3

        t, p, dpp, ipi, ipe, _, _, q, dpq, iqi, iqe, _, _, _ = airburst(
            groundRange, height, Y, None, False
        )

        P[j][i] = p / 1e6
        Q[j][i] = q / 1e6

        T[j][i] = t

        if ipi < 1:
            ipi = 1

        if iqi < 1:
            iqi = 1

        IP[j][i] = ipi
        DPP[j][i] = dpp
        IQ[j][i] = iqi
        DPQ[j][i] = dpq

        j += 1

    i += 1

IP = ma.masked_where(IP is None, IP)
IQ = ma.masked_where(IQ is None, IQ)
DPP = ma.masked_where(DPP is None, DPP)
DPQ = ma.masked_where(DPQ is None, DPQ)
IP = IP / 1e3
IQ = IQ / 1e3


fig, axs = plt.subplots(2, 3)
axu, axl = axs
ax4, ax5, ax6 = axu
ax1, ax2, ax3 = axl
highs = (ax4, ax5, ax6)
lows = (ax1, ax2, ax3)


mode = 0  # mode 0: overpressure mode 1: dynamic pressure


def onclick(event):
    global mode

    if mode == 0:
        cmap = "Spectral_r"
    else:
        cmap = "coolwarm"
    plt.clf()

    axs = fig.subplots(2, 3)

    axu, axl = axs
    ax4, ax5, ax6 = axu
    ax1, ax2, ax3 = axl
    highs = (ax4, ax5, ax6)
    lows = (ax1, ax2, ax3)

    for ax in lows:
        ax.set_xlim(0, xmed / 1e3)
        ax.set_ylim(0, ymed / 1e3)
        ticks = tuple(i * imed / 1e3 for i in range(int(xmed // imed) + 1))
        ax.xaxis.set_ticks(ticks)
        ax.yaxis.set_ticks(ticks)
        ax.grid(which="major", linestyle=":", color="#DDDDDD", linewidth=1)

    for ax in highs:
        ax.set_xlim(0, xmax / 1e3)
        ax.set_ylim(0, ymax / 1e3)
        ax.xaxis.set_ticks(tuple(i * imax / 1e3 for i in range(int(xmax // imax) + 1)))
        ax.grid(which="major", linestyle="-", color="#DDDDDD", linewidth=1)

    ax2.yaxis.tick_right()
    ax2.yaxis.set_label_position("right")
    ax3.yaxis.tick_right()
    ax3.yaxis.set_label_position("right")

    ax5.yaxis.tick_right()
    ax5.yaxis.set_label_position("right")
    ax6.yaxis.tick_right()
    ax6.yaxis.set_label_position("right")

    fig.suptitle("{:.0f}kT Airburst".format(Y), fontsize="x-large")
    ax4.set_title("Overpressure Damaging Effects")
    if mode == 0:
        ax5.set_title("Overpressure Arrival & Positive Phase Duration (s)")
        ax6.set_title("Overpressure Impulse (kPa-s)")
    else:
        ax5.set_title("Dynamic Pressure Arrival & Positive Phase Duration (s)")
        ax6.set_title("Dynamic Pressure Impulse (kPa-s)")

    for ax in lows + highs:
        ax.set_aspect("equal")
        ax.set_xlabel("Ground Range (km)")
        ax.set_ylabel("Burst Height (km)")

    for ax in highs:
        rect = patches.Rectangle(
            (0, 0),
            xmed / 1e3,
            ymed / 1e3,
            linewidth=1.5,
            edgecolor="r",
            facecolor="none",
        )
        ax.add_patch(rect)

    dynamic = (0.12, 0.14, 0.50, 0.70)
    dynamic = tuple(i * 0.098 for i in dynamic)
    csq = ax1.contour(
        x,
        y,
        Q,
        levels=dynamic,
        colors="darkorange",
        linestyles=[
            (0, (5, 1)),
            (0, (5, 1)),
            "solid",
            "solid",
        ],
    )

    personnel = (0.17, 0.20, 0.30, 0.35, 0.53, 0.55, 0.95, 1)
    personnel = tuple(i * 0.098 for i in personnel)
    csp = ax1.contour(
        x,
        y,
        P,
        levels=personnel,
        colors="black",
        linestyles=(
            "dotted",
            "dotted",
            "dashed",
            "dashed",
            (0, (5, 1)),
            (0, (5, 1)),
            "solid",
            "solid",
        ),
    )

    nmp, lblp = csp.legend_elements()
    nmq, lblq = csq.legend_elements()
    lblp = ("Light", "Medium", "Heavy", "Extreme")
    lblq = ("Hvy. (Dyn.)", "Ex. (Dyn.)")

    csph = ax1.contourf(x, y, P, levels=personnel, colors=("lightgrey", "none"))

    i = 0
    for h in (1004, 880, 676, 595, 470, 445, 335, 310):
        ax1.text(
            1 * Y3 / 1e3,
            h * Y3 / 1e3,
            "{:.0f}kPa".format(personnel[i] * 1e3),
            va=("bottom", "top")[i % 2],
            ha="left",
            alpha=0.75,
        )
        i += 1

    legend1 = ax1.legend(
        nmp[::2] + nmq[::2],
        lblp + lblq,
        title="Exposed Personnel",
        ncol=1,
        loc="upper right",
    )

    shelter = (0.04, 0.18, 0.40)
    shelter = tuple(i * 0.098 for i in shelter)

    css = ax4.contour(
        x,
        y,
        P,
        levels=shelter,
        colors="seagreen",
        linestyles=("dotted", "dashed", "solid"),
    )
    nms, lbls = css.legend_elements()
    lbls = ("Light", "Medium", "Severe")
    legend2 = ax4.legend(nms, lbls, title="Low Buildings", ncol=1, loc="upper right")

    i = 0
    for h in (3025, 950, 545):
        ax4.text(
            1 * Y3 / 1e3,
            h * Y3 / 1e3,
            "{:.0f}kPa".format(shelter[i] * 1e3),
            va="bottom",
            color="seagreen",
            ha="left",
            alpha=0.75,
        )
        i += 1

    def t(y, x):
        t, _, _, _, _, _, _, _, _, _, _, _, _, _ = airburst(x, y, Y, None, False)
        return t

    cst2 = ax2.contour(
        x,
        y,
        T,
        colors="black",
        levels=[
            t(0, x * imed)
            for x in range(1, int(sqrt(xmed**2 + ymed**2) // imed + 2))
        ],
        extend="min",
        vmin=0,
        vmax=t(ymax, xmax),
        # linestyles=":",
    )
    ax2.clabel(cst2, fontsize="smaller", fmt=lambda x: "{:.2f} s".format(x))

    cst5 = ax5.contour(
        x,
        y,
        T,
        colors="black",
        levels=[
            t(0, x * imax)
            for x in range(1, int(sqrt(xmax**2 + ymax**2) // imax) + 2)
        ],
        extend="both",
    )
    ax5.clabel(cst5, fontsize="smaller", fmt=lambda x: "{:.1f} s".format(x))

    if mode == 0:
        pd = DPP
    else:
        pd = DPQ

    pdmin = np.nanmin(pd)
    pdmax = np.nanmax(pd)

    csd2 = ax2.contourf(
        x,
        y,
        pd,
        cmap=cmap,
        alpha=0.75,
        levels=np.linspace(pdmin, pdmax, 11),
    )

    _, top = ax2.get_ylim()
    axins2 = inset_axes(
        ax2,
        width="100%",
        height="100%",
        loc="lower right",
        borderpad=0,
        bbox_to_anchor=(-1.25 * top / 10, 0, top / 10, top),
        bbox_transform=ax2.transData,
    )

    cbar2 = fig.colorbar(
        csd2, cax=axins2, format="%.2f", label="Positive Phase Duration (s)"
    )

    csd5 = ax5.contourf(
        x,
        y,
        pd,
        cmap=cmap,
        alpha=0.75,
        linestyles=":",
        levels=np.linspace(pdmin, pdmax, 11),
    )
    _, top = ax5.get_ylim()
    axins5 = inset_axes(
        ax5,
        width="100%",
        height="100%",
        loc="lower right",
        borderpad=0,
        bbox_to_anchor=(-1.25 * top / 10, 0, top / 10, top),
        bbox_transform=ax5.transData,
    )
    cbar5 = fig.colorbar(
        csd5, cax=axins5, format="%.2f", label="Positive Phase Duration (s)"
    )

    if mode == 0:
        ip = IP
    elif mode == 1:
        ip = IQ

    ipmin = np.nanmin(ip)
    ipmax = np.nanmax(ip)

    csd3 = ax3.contourf(
        x,
        y,
        ip,
        alpha=0.75,
        cmap=cmap,
        locator=ticker.LogLocator(),
        levels=np.logspace(log10(ipmin), log10(ipmax), 11),
        # levels=np.linspace(ipmin, ipmax, 10),
    )

    _, top = ax3.get_ylim()
    axins3 = inset_axes(
        ax3,
        width="100%",
        height="100%",
        loc="lower left",
        borderpad=0,
        bbox_to_anchor=(-1.25 * top / 10, 0, top / 10, top),
        bbox_transform=ax3.transData,
    )
    cbar3 = fig.colorbar(csd3, cax=axins3, format="%.2f", label="Impulse (kPa-s)")

    csd6 = ax6.contourf(
        x,
        y,
        ip,
        alpha=0.75,
        cmap=cmap,
        locator=ticker.LogLocator(),
        levels=np.logspace(log10(ipmin), log10(ipmax), 11),
    )

    _, top = ax6.get_ylim()
    axins6 = inset_axes(
        ax6,
        width="100%",
        height="100%",
        loc="lower left",
        borderpad=0,
        bbox_to_anchor=(-1.25 * top / 10, 0, top / 10, top),
        bbox_transform=ax6.transData,
    )

    cbar6 = fig.colorbar(csd6, cax=axins6, format="%.2f", label="Impulse (kPa-s)")

    axins2.yaxis.set_ticks_position("left")
    axins2.yaxis.set_label_position("left")
    axins3.yaxis.set_ticks_position("left")
    axins3.yaxis.set_label_position("left")
    axins5.yaxis.set_ticks_position("left")
    axins5.yaxis.set_label_position("left")
    axins6.yaxis.set_ticks_position("left")
    axins6.yaxis.set_label_position("left")
    cbar2.ax.tick_params(labelsize="small")
    cbar3.ax.tick_params(labelsize="small")
    cbar5.ax.tick_params(labelsize="small")
    cbar6.ax.tick_params(labelsize="small")

    plt.draw()  # redraw

    mode = (mode + 1) % 2


if __name__ == "__main__":
    fig.canvas.mpl_connect("button_press_event", onclick)
    onclick(None)
    plt.show()
