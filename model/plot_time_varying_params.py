import matplotlib.pyplot as plt
import numpy as np

from deampy.parameters import Constant, TimeDependentSigmoid
from deampy.plots.plot_support import output_figure
from definitions import ROOT_DIR

SIM_DURATION = 25
T0 = 0
X_RANGE = (-0.1, 26)
Y_RANGE = (0, 1.05)


def sigmoid(t, b, f_min, f_max, t_mid, t_min=0):
    # sigmoid function
    par = TimeDependentSigmoid(
        par_b=Constant(b),
        par_t_min=Constant(t_min),
        par_t_middle=Constant(t_mid),
        par_min=Constant(f_min),
        par_max=Constant(f_max)
    )
    val = par.sample(time=t)
    return val if val > 0 else None


def plot_sigmoid_functions(b, f_min, f_max, t_mid, t_min,
                           bs, f_mins, f_maxs, t_mids, t_mins,
                           x_min_max, y_label, x_label,
                           vertical_line=None, fig_size=(7.5, 3.2), round_b=None, file_name='sigmoid.png'):

    # ------------------
    ts = np.linspace(start=x_min_max[0], stop=x_min_max[1])
    fs = []
    legends = []
    titles = []

    # varying b
    if bs is not None:
        titles.append('Varying ' + r'$b$')
        legs = []  # legends
        ys = []
        for v in bs:
            if round_b is None:
                legs.append(r'$b=$' + str(v))
            else:
                legs.append(r'$b=$' + str(round(v, round_b)))
            ys.append([sigmoid(t, b=v, t_mid=t_mid, f_min=f_min, f_max=f_max, t_min=t_min) for t in ts])
        legends.append(legs)
        fs.append(ys)

    # varying b_min
    if f_mins is not None:
        titles.append('Varying ' + r'$b_{min}$')
        legs = []  # legends
        ys = []
        for v in f_mins:
            legs.append(r'$b_{min}=$' + str(v))
            ys.append([sigmoid(t, b=b, t_mid=t_mid, f_min=v, f_max=f_max, t_min=t_min) for t in ts])
        legends.append(legs)
        fs.append(ys)

    # varying b_max
    if f_maxs is not None:
        titles.append('Varying ' + r'$b_{max}$')
        legs = []  # legends
        ys = []
        for v in f_maxs:
            legs.append(r'$b_{max}=$' + str(v))
            ys.append([sigmoid(t, b=b, t_mid=t_mid, f_min=f_min, f_max=v, t_min=t_min) for t in ts])
        legends.append(legs)
        fs.append(ys)

    # varying t_mid
    if t_mids is not None:
        titles.append('Varying ' + r'$t_{mid}$')
        legs = []  # legends
        ys = []
        for v in t_mids:
            legs.append(r'$t_{mid}=$' + str(v))
            ys.append([sigmoid(t, b=b, t_mid=v, f_min=f_min, f_max=f_max, t_min=t_min) for t in ts])
        legends.append(legs)
        fs.append(ys)

    # varying t_mins
    if t_mins is not None:
        titles.append('Varying ' + r'$t_{min}$')
        legs = []  # legends
        ys = []
        for v in t_mins:
            legs.append(r'$t_{min}=$' + str(v))
            ys.append([sigmoid(t, b=b, t_mid=t_mid, f_min=f_min, f_max=f_max, t_min=v) for t in ts])
        legends.append(legs)
        fs.append(ys)

    # plot
    fig, axarr = plt.subplots(1, len(fs), sharey=True, figsize=fig_size)
    for i, ax in enumerate(axarr):
        ax.set_title(titles[i])

        for j in range(3):
            ax.plot(ts, fs[i][j], label=legends[i][j])  # color='b', linestyle='-')

        ax.set_ylim(Y_RANGE)
        ax.set_xlim(X_RANGE)
        if vertical_line is not None:
            ax.axvline(x=vertical_line, c='k', linestyle='--', linewidth=0.5)
        ax.set_xlabel(x_label)
        ax.legend(fontsize='x-small') # loc=2

    axarr[0].set_ylabel(y_label)
    plt.tight_layout()
    output_figure(plt, filename=file_name)
    plt.show()


# probability of novel strain over time
plot_sigmoid_functions(b=0.3, bs=[0.1, 0.3, 0.5],
                       f_min=0.8, f_mins=[0.7, 0.8, 0.9],
                       f_max=1, f_maxs=None,
                       t_mid=10, t_mids=[7, 10, 13],
                       t_min=0, t_mins=None,
                       y_label=r'$\gamma(t)$',
                       x_min_max=[T0, SIM_DURATION],
                       x_label='Simulation year {}'.format(r'$(t)$'),
                       fig_size=(6, 2.8),
                       file_name=ROOT_DIR+'/analysis/figures/params/fitness.png')
