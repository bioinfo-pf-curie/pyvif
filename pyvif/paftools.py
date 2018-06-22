# -*- coding: utf-8 -*-

"""PAF toolbox"""

import matplotlib.pyplot as plt
import pandas as pd

from . import logger, _PAF_COLNAMES
from .bamtools import bam_to_paf
from .exception import BadInputException
from .plots import plot_histogram


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
            raise BadInputException(msg)
        except AttributeError:
            pass
        # Check bam file
        try:
            return bam_to_paf(filin)
        except ValueError:
            pass
        msg = "No correct input provided."
        raise BadInputException(msg)

    def number_mapped_reads(self):
        """ Return the number of mapped reads.
        """
        return len(self.df['q_name'].unique())

    def plot_length(self, targets=None, filename=None, bins=100):
        """ Mapped reads length histogram.

        :params list targets: set of targets contigs name.
        :params str filename: filename to save images.
        :params int bins: number of bins for barplot.

        Return image filename or show the plot.
        """
        # set plot title and get reads of interest
        title = "Mapped read lengths"
        if targets:
            title += " on {} histogram".format(", ".join(targets))
            paf = self.df.loc[self.df['chr'].isin(targets)]
        else:
            title += " histogram"
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
            bins=bins
        )

    def plot_positions(self, filename=None):
        """Barplot of reads position in reference.

        :params str filename: filename to save images.

        Return image filename or dictionnary with barplot results.
        """
        # generate data
        position = self.df.groupby('chr').q_name.nunique()
        # create figure
        plt.figure(figsize=(10, 7.5))
        ax = plt.subplot(111)

        # remove top and right frame lines
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

        # increase font size and lean x labels
        plt.xticks(fontsize=14, rotation=50)
        plt.yticks(fontsize=14)
        plt.xlabel("Contigs", fontsize=16)
        plt.ylabel("Count", fontsize=16)
        for tick in ax.xaxis.get_majorticklabels():
            tick.set_horizontalalignment('right')
        plt.title("Positions of mapped reads", fontsize=18)

        plt.bar(position.index, position, color="#3F5D7D")

        if filename:
            plt.savefig(filename, bbox_inches="tight")
            return filename
        plt.show()

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
        plt.figure(figsize=(10, 7.5))
        ax = plt.subplot(111)

        # increase font size
        plt.xticks(fontsize=14)
        plt.yticks(fontsize=14)
        plt.xlabel("Read lengths", fontsize=16)
        plt.ylabel("Percentage aligned", fontsize=16)

        paf = self.df.copy()
        paf['a_length'] = paf['r_end'] - paf['r_start']
        paf['%_aligned'] = (paf['a_length'] / paf['q_length']) * 100
        if targets:
            paf = paf.loc[paf['chr'].isin(targets)]

        plt.hist2d(paf['q_length'], paf['%_aligned'], bins=bins,
                   cmap=cm.afmhot_r)
        plt.colorbar()
        if filename:
            plt.savefig(filename, bbox_inches="tight")
            return filename
        plt.show()

    def plot_number_pass(self, filename=None):
        """ Plot the distribution of pass number.

        :params str filename: filename to save image.

        Return image filename or show the plot in ipython.
        """

        paf = self.df.copy()
        name, zmw, interval = paf['q_name'].str.split('/').str
        paf['zmw'] = pd.to_numeric(zmw, errors='coerce')
        # count number of pass for each zmw
        zmw_count = paf.groupby('q_name').zmw.max().value_counts()
        # count number of number of pass
        pass_count = zmw_count.value_counts()

        # create figure
        plt.figure(figsize=(10, 7.5))
        ax = plt.subplot(111)

        # remove top and right frame lines
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

        # increase font size and lean x labels
        plt.xticks(fontsize=14)
        plt.yticks(fontsize=14)
        plt.xlabel("Number of pass", fontsize=16)
        plt.ylabel("Count", fontsize=16)
        plt.title("Number of number of pass", fontsize=18)

        plt.bar(pass_count.index, pass_count, color="#3F5D7D")

        if filename:
            plt.savefig(filename, bbox_inches="tight")
            return filename
        plt.show()
