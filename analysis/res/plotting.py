"""Histogram plotting module for Oxford HH->4B pheno code

Import this module (``from plotting import *``) then use the `plot1D`
`plot2D` functions. Loading Python and Matplotlib takes a while, so
try to do all your plotting in one Python script.

Functions:
plot1D -- Plot 1D histograms. Allows you to plot multiple 1D histograms
          on the same axes. All histograms are plotted normalized by area.
plot2D -- Plot a single 2D histogram as a heatmap.
"""
from __future__ import division, print_function
import matplotlib.pyplot as plt
import numpy as np
import numexpr as ne
import seaborn as sns
from itertools import cycle

__all__ = ["plot1D", "plot2D"]

plt.switch_backend("PDF")


def plot1D(histograms, xlabel='Observable bin', title='', output=''):
    """Plot 1D histograms

    Arguments:
    histograms -- Dict with keys being legend labels (use mathtext) and
                  values being histogram data filenames
    xlabel -- X axis label
    title  -- Plot title
    output -- Output filename (file will be saved under `plots` directory)
    """
    col_cycle = cycle(['C0', 'C1', 'C2', 'C3'])
    fig = plt.figure(figsize=(8, 8))
    ax = plt.gca()
    ax.set_xlabel(xlabel)
    ax.set_title(title)
    ax.set_ylabel("Arbitrary units")
    for (label, filename), color in zip(histograms.items(), col_cycle):
        # Only read out errminus, because the two errors are the same anyway
        xlow, xhigh, val, err, _ = np.loadtxt(
            filename, skiprows=4, unpack=True)
        errLow = val - err
        errHigh = val + err
        # normalization
        normalization = ne.evaluate("sum(val / (xhigh - xlow))")
        x = np.stack((xlow, xhigh)).ravel(1)  # interleaves xlow and xhigh
        y = np.stack((val, val)).ravel(1)  # doubles up val
        y /= normalization
        yerrLow = np.stack((errLow, errLow)).ravel(1)
        yerrHigh = np.stack((errHigh, errHigh)).ravel(1)
        yerrLow /= normalization
        yerrHigh /= normalization

        ax.fill_between(x, yerrLow, yerrHigh, color=color, alpha=0.5)
        ax.plot(x, y, label=label, color=color)

    ax.xaxis.grid(True)
    ax.yaxis.grid(True)

    leg = ax.legend(loc='best')
    leg.get_frame().set_alpha(0.8)
    fig.savefig('plots/' + output)
    plt.close(fig)


def plot2D(title, filename, xlabel, ylabel, output):
    """Plot a 2D histogram

    Arguments:
    title -- Plot title
    filename -- Histogram data filename
    xlabel -- X axis label
    ylabel -- Y axis label
    output -- Output filename (file will be saved under `plots` directory)
    """
    fig = plt.figure(figsize=(8, 8))
    ax = plt.gca()
    xlow, xhigh, ylow, yhigh, val, _, _ = np.loadtxt(
        filename, skiprows=4, unpack=True)
    length = xlow.size
    xlowNew = np.unique(xlow)
    xhighNew = np.unique(xhigh)
    if (xlowNew.size != xhighNew.size):
        print("PLOT2D ERR: xlowNew.size != xhighNew.size")
        return
    xlength = xlowNew.size
    if (length % xlength != 0):
        print("PLOT2D ERR: y bins non-constant")
    # Above checks are just-in-case. Shouldn't be necessary
    ylength = length // xlength
    ylowNew = ylow[0:ylength]  # noqa
    yhighNew = yhigh[0:ylength]  # noqa

    x = ne.evaluate("xlowNew + (xhighNew - xlowNew) / 2")
    y = ne.evaluate("ylowNew + (yhighNew - ylowNew) / 2")
    z = val.reshape((xlength, ylength)).T
    sns.heatmap(
        z, xticklabels=x, yticklabels=y, ax=ax, square=True, linewidths=0.5)
    plt.xticks(rotation=90)
    plt.yticks(rotation=0)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    fig.savefig('plots/' + output)
    plt.close(fig)
