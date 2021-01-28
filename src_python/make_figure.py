import math
import numpy as np
import matplotlib.pyplot as plt
import mpl_toolkits.axisartist as axisartist


def draw():
    fig = plt.figure(dpi=192)
    ax = axisartist.Subplot(fig, 111)
    fig.add_axes(ax)

    # 上軸と右軸を表示しない
    ax.axis["right"].set_visible(False)
    ax.axis["top"].set_visible(False)

    delta = 0.025
    xrange = np.arange(-2, 2, delta)
    yrange = np.arange(-2, 2, delta)
    X, Y = np.meshgrid(xrange, yrange)

    for i in range(5):
        rad = math.pi*(i+1)*15/180
        x1 = math.cos(rad)
        y1 = math.sin(rad)
        ax.annotate('', xy=[x1*0.5, y1*0.5], xytext=[x1, y1],
                    arrowprops=dict(shrink=0, width=1, headwidth=8,
                                    headlength=10, connectionstyle='arc3',
                                    facecolor='red', edgecolor='red')
                    )
        x2 = 1.5*math.cos(rad)
        y2 = 1.5*math.sin(rad)
        ax.plot([x1, x2], [y1, y2], color="blue")

    points = [[], []]
    for i in range(45):
        rad = math.pi*i*2/180
        x = math.cos(rad)
        y = math.sin(rad)
        points[0].append(x)
        points[1].append(y)
        x *= 1.25
        y *= 1.25
        points[0].append(x)
        points[1].append(y)
        x *= 1.2
        y *= 1.2
        points[0].append(x)
        points[1].append(y)
    ax.plot(points[0], points[1], marker='o', markersize=2, color="black", linestyle='None')

    # 軸の設定
    ax.set_xlim(0, 2)
    ax.set_ylim(0, 2)

    # アスペクト比の設定
    ax.set_aspect('equal', adjustable='box')

    # ax.tick_params(labelbottom=False,
    #                labelleft=False,
    #                labelright=False,
    #                labeltop=False)

    ax.tick_params(bottom=False, left=False, right=False, top=False)
    for direction in ["left",  "bottom"]:
        ax.axis[direction].toggle(all=False)
        ax.axis[direction].set_axisline_style("-|>")

    fig.savefig("test.png")


if __name__ == "__main__":
    draw()
