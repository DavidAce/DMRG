import itertools
import seaborn as sns
from matplotlib import rc

rc('text', usetex=True)
sns.set_style(rc={'text.usetex': True})


# sns.set_palette(sns.color_palette("colorblind", 8))
# ed_palette = sns.color_palette("bright", 3)
# dmrg_palette = sns.color_palette("colorblind", 5)
# sns.set_palette(sns.color_palette("muted", 8))
# sns.set_palette(sns.husl_palette(4, h=0.6, s=0.9, l=0.5))
# sns.set_palette(sns.color_palette("Accent", n_colors=10))
# sns.set_palette(sns.choose_colorbrewer_palette("qualitative"))

####################################################################################
####### Use these settings for presentations or where big fonts are needed #########
####################################################################################
# rc('font', **{'family': 'serif', 'serif': ['Palatino']})
# rc('font', **{'family': 'sans-serif', 'sans-serif': ['Helvetica']})
# rc('text', usetex=True)
# # sns.set(style="whitegrid", font_scale=1.5 ,rc={'text.usetex' : True})
# sns.set(style="darkgrid", font_scale=1.5, font='Helvetica',rc={"lines.linewidth": 1.2})
# paper_rc = {'lines.linewidth': 2, 'lines.markersize': 16}
# sns.set_style({"axes.facecolor": ".9"},rc={'text.usetex' : True})
# sns.set_palette(sns.color_palette("husl", 5))
####################################################################################


def get_markerstyles():
    return itertools.cycle(('o', 'v', '^', '<', '>', 'p', '8', '*', 'h', '+', 'x', 'D'))


def get_linestyles():
    return itertools.cycle(('-', '--', '-.', ':'))


def plt_scatter(x, y, xlabel, ylabel, ax, markerlist, label='', xscale='linear', yscale='linear', markersize=20):
    ax.scatter(x, y, marker=next(markerlist), label=label, alpha=.90, s=markersize)
    if xlabel != '':  ax.set_xlabel(xlabel)
    if ylabel != '':  ax.set_ylabel(ylabel)
    ax.set_xscale(xscale)
    ax.set_yscale(yscale)
    if label != '': ax.legend()
