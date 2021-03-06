# -*- coding: utf-8 -*-

"""PAF toolbox"""

import matplotlib.pyplot as plt
import pandas as pd

from pyvif import logger, _PAF_COLNAMES
from pyvif.bamtools import bam_to_paf
from pyvif.plots import plot_histogram, init_plot, generate_plot
from pyvif.config import COLOR


class PAF(object):
    """ Pairwise mApping Format.
    PAF is a text format describing the approximate mapping positions between
    two set of sequences.
        - chr: chromosome name
        - r_start: start position on reference
        - r_end: end position on reference
        - mapq: mapping quality score
        - strand: sens of mapping
        - q_name: read name
        - q_length: read length
        - q_start: strat position on read
        - q_end: end position on read
    """
    def __init__(self, paf):
        """.. rubric: constructor

        :param paf: pd.DataFrame in paf format returned by
            :meth:`bamtools.bam_to_paf` or a bam file.
        """
        self.df = self._get_df(paf)

    def _get_df(self, filin):
        # Check DataFrame
        try:
            if filin.columns.tolist() == _PAF_COLNAMES:
                return filin
            msg = ("DataFrame input must have these columns in this order:"
                   " {}".format(" - ".join(_PAF_COLNAMES)))
            logger.error(msg)
            raise TypeError(msg)
        except AttributeError:
            pass
        # Check bam file
        try:
            return bam_to_paf(filin)
        except ValueError:
            pass
        msg = "No correct input provided."
        raise TypeError(msg)

    def reads_count(self):
        """ Return the number of mapped reads.
        """
        return len(self.df['q_name'].unique())

    def plot_length(self, filename=None, targets=None, title=None, bins=100):
        """ Mapped reads length histogram.

        :params list targets: set of targets contigs name.
        :params str filename: filename to save images.
        :params int bins: number of bins for barplot.

        Return image filename or show the plot.
        """
        # set plot title and get reads of interest
        if title is None:
            title = "Mapped read lengths"
            if targets:
                title += " on {} histogram".format(", ".join(targets))
            else:
                title += " histogram"

        # select targets
        try:
            paf = self.df.loc[self.df['chr'].isin(targets)]
        except TypeError:
            paf = self.df

        # get data
        lengths = paf.groupby("q_name").q_length.max()
        
        # plot histogram
        plot_histogram(
            data=lengths,
            title=title,
            xlabel="Read lengths",
            ylabel="Count",
            filename=filename,
            bins=bins,
            log=True
        )

    def plot_positions(self, filename=None):
        """Barplot of reads position in reference.

        :params str filename: filename to save images.

        Return image filename.
        """
        # generate data
        position = self.df.groupby('chr').q_name.nunique()

        # init plot
        fig, ax = init_plot("Positions of mapped reads", 'Contigs', 'Count')
        ax.set_yscale('log')
        plt.xticks(fontsize=14, rotation=50)
        for tick in ax.xaxis.get_majorticklabels():
            tick.set_horizontalalignment('right')
        plt.bar(position.index, position, color=COLOR)
        return generate_plot(filename, fig)

    def plot_length_vs_mapping(self, targets=None, filename=None, bins=100):
        """ 2D histogram of target positive reads length versus percentage of
        reads that maps on reference.

        :params set targets: set of targets contigs name.
        :params str filename: filename to save images.
        :params int Nlevel: number of level for colors.
        :params int bins: number of bins for barplot.

        Return image filename or show the plot in ipython.
        """
        from matplotlib import cm
        fig, ax = init_plot("", "Read lengths", "Percentage aligned")

        paf = self.df.copy()
        paf['a_length'] = paf['r_end'] - paf['r_start']
        paf['%_aligned'] = (paf['a_length'] / paf['q_length']) * 100
        if targets:
            paf = paf.loc[paf['chr'].isin(targets)]

        plt.hist2d(paf['q_length'], paf['%_aligned'], bins=bins,
                   cmap=cm.viridis)
        plt.colorbar()
        return generate_plot(filename, fig)


    def plot_number_pass(self, filename=None):
        """ Plot the distribution of pass number.

        :params str filename: filename to save image.

        Return image filename or show the plot in ipython.
        """

        paf = self.df.copy()
        try:
            name, zmw, interval = paf['q_name'].str.split('/').str
        except ValueError:
            logger.warning("You did not use subreads."
                           " This plot is unavailable.")
            return None
        paf['zmw'] = pd.to_numeric(zmw, errors='coerce')
        # count number of pass for each zmw
        zmw_count = paf.groupby('q_name').zmw.max().value_counts()
        # count number of number of pass
        pass_count = zmw_count.value_counts()

        fig, ax = init_plot("Count number of pass", "Number of pass", 'Count')
        ax.set_yscale('log')
        plt.bar(pass_count.index, pass_count, color=COLOR)
        return generate_plot(filename, fig)
