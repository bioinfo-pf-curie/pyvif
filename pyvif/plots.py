# -*- coding: utf-8 -*-

"""Plot toolbox"""

import matplotlib.pyplot as plt


def plot_histogram(data, title, xlabel, ylabel, filename=None, bins=100,
                   color='#3F5D7D'):
    """ Plot histogram.
    """
    plt.figure(figsize=(10, 7.5))
    ax = plt.subplot(111)

    # remove top and right frame lines
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()

    # increase label font size
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.xlabel(xlabel, fontsize=16)
    plt.ylabel(ylabel, fontsize=16)
    plt.title(title, fontsize=18)

    plt.hist(data, bins=bins, color=color)
    if filename:
        plt.savefig(filename, bbox_inches="tight")
        return filename
    plt.show()