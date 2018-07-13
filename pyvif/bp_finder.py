# -*- coding: utf-8 -*-

""" A second version of breakpoint detection

I detect breakpoint with the minimap2 mapping and breakpoints sequences.
"""

from collections import OrderedDict
from itertools import takewhile

import pandas as pd
import matplotlib.pyplot as plt

from . import logger
from .paftools import PAF
from .plots import init_plot

_TRANSTAB = str.maketrans("ACGTacgt", "TGCAtgca")


class BreakpointFinder(object):
    """ Object that detect all breakpoint present in PacBio data.
    It uses a mapping on human and viruses genomes.
    """
    def __init__(self, human, virus):
        """.. rubric:: constructor

        :param human: PAF* format of mapping on human genome or BAM file.
        :param virus: PAF* format of mapping on virus genomes or BAM file.

        *PAF correspond to :class:`pbcapture.bamtools.PAF` object or a PAF
        file created by :meth:`pbcapture.bamtools.PAF.to_csv`.
        """
        self.paf, self.virus_contigs = self._init_paf(human, virus)
        self.palindromics = []
        self.bps = self.find_breakpoints()

    def _init_paf(self, human, virus):
        """ Initiate intersect paf with alignment on human and on target.
        Only reads that maps on virus genomes are kept.

        Return merged PAFs as a :class:`pd.DataFrame` and contigs of virus.
        """
        human_paf = self._create_paf_format(human)
        virus_paf = self._create_paf_format(virus)

        # get human alignment of reads HPV+
        virus_reads_name = virus_paf['q_name'].unique()
        human_align = human_paf.loc[human_paf['q_name'].isin(virus_reads_name)]

        paf = pd.concat([human_align, virus_paf])
        paf = paf.sort_values(["q_name", "q_start"])
        paf = paf.set_index("q_name")
        return paf, virus_paf["chr"].unique().tolist()

    def _create_paf_format(self, filin):
        """ The class can be init with BAM file or PAF format generated by
            :class:`pbcapture.bamtools.PAF`.

            Return the PAF as a :class:`pd.DataFrame`.
        """
        # Check DataFrame and bamfile
        try:
            return PAF(filin).df
        except TypeError:
            pass
        # Check PAF format
        try:
            filin.reads_count()
        except AttributeError:
            pass
        msg = "No correct input provided."
        raise TypeError(msg)

    def find_breakpoints(self):
        """ Find all breakpoints between the reference and the contig of
        virus.
        """
        logger.info("Breakpoint research is running...")
        self.palindromics = []
        breakpoints = list()
        for read_name in self.paf.index.unique():
            try:
                breakpoints += self._get_read_breakpoints(read_name)
            except TypeError:
                pass
        logger.info("{} breakpoints are found.".format(len(breakpoints)))
        return pd.DataFrame(breakpoints)

    def _get_read_breakpoints(self, name):
        try:
            read_paf = self.paf.loc[name]
        except KeyError:
            return None
        # zip function is faster than itertuples method
        # try if read_paf is a dataframe and not a Series
        # If it is a Serie -> read have only one alignment
        try:
            iter_align = zip(
                read_paf["chr"], read_paf["r_start"], read_paf["r_end"],
                read_paf["q_start"], read_paf["q_end"], read_paf["strand"],
                read_paf["mapq"]
            )
        except TypeError:
            return None

        # init previous variables
        prev_chr = prev_coord = prev_strand = prev_mapq = None
        bp_list = list()

        # find breakpoint in read
        for cur_chr, *cur_coord, cur_strand, cur_mapq in iter_align:
            if cur_chr in self.virus_contigs:
                if prev_chr != cur_chr:
                    try:
                        bp_list.append(
                            self._get_breakpoint(
                                prev_chr, prev_coord, prev_strand,
                                prev_mapq, cur_chr, cur_coord,
                                cur_strand, name
                            )
                        )
                    except TypeError:
                        pass
            elif prev_chr in self.virus_contigs:
                try:
                    bp_list.append(
                        self._get_breakpoint(
                            cur_chr, cur_coord, cur_strand, cur_mapq, prev_chr,
                            prev_coord, prev_strand, name
                        )
                    )
                except TypeError:
                    pass
            (prev_chr, prev_coord, prev_strand, prev_mapq) = (
                cur_chr, cur_coord, cur_strand, cur_mapq)

        # Check if the read is palindromics
        if _read_is_palindromic(bp_list):
            self.palindromics.append(bp_list)
            return []
        return bp_list

    def _get_breakpoint(self, human_chr, human_coord, human_strand, human_mapq,
                        virus_chr, virus_coord, virus_strand, read_name):
        """ Create breakpoint dictionnary to create dataframe.
        """
        # check if hg - hpv or hpv - hg case
        i = 0 if human_coord[2] < virus_coord[2] else 1

        # get coord if strand is '+'
        too_far = abs(human_coord[3 - i] - virus_coord[2 + i]) > 20
        if human_mapq > 20 and not too_far:
            bpstart_human = human_coord[1 - i]
            end_human = human_coord[0 + i]
        else:
            human_strand = human_chr = bpstart_human = end_human = "Unknown"
        bpstart_virus = virus_coord[0 + i]
        end_virus = virus_coord[1 - i]

        # switch start / end if strand is '-'
        if human_strand == '-':
            bpstart_human, end_human = end_human, bpstart_human
        if virus_strand == '-':
            bpstart_virus, end_virus = end_virus, bpstart_virus

        if virus_strand == '-':
            # the virus strand always set as '+'
            virus_strand = '+'
            if human_strand != 'Unknown':
                human_strand = '+' if human_strand == '-' else '-'

        return OrderedDict([
            ("chromosome", human_chr),
            ("bpstart_human", bpstart_human),
            ("end_human", end_human),
            ("strand_human", human_strand),
            ("virus_contig", virus_chr),
            ("bpstart_virus", bpstart_virus),
            ("end_virus", end_virus),
            ("strand_virus", virus_strand),
            ("read", read_name),
        ])

    def clustering_breakpoints(self, human_thd=25, virus_thd=25):
        """ Generate a clustering with bpstart_human and with bpstart_virus.

        :params int threshold: maximal distance between two breakpoint start.

        Add 'human_clust' and 'virus_clust' in dataframe.
        """
        # clustering with human breakpoint start and chromosome
        known_bp = self.bps.loc[self.bps['chromosome'] != 'Unknown']
        clustering = known_bp.sort_values('bpstart_human')\
                             .groupby('chromosome')['bpstart_human']\
                             .diff().gt(human_thd).cumsum()
        self.bps.loc[clustering.index, 'human_clust'] = clustering
        grouped = self.bps.dropna().groupby(['chromosome', 'human_clust'])
        logger.info("They are {} different breakpoints in human genome."
                    .format(len(grouped.groups)))
        # clustering with virus breakpoint start
        clustering = self.bps.sort_values('bpstart_virus')['bpstart_virus']\
                             .diff().gt(virus_thd).cumsum()
        self.bps.loc[clustering.index, 'virus_clust'] = clustering
        logger.info("They are {} different breakpoints in virus genome."
                    .format(len(self.bps['virus_clust'].unique())))

    def summarize_human_clustering(self, threshold=20):
        """ Return a dataframe that summarizes bigger clusters.

        :params int threshold: minimum size of cluster.
        """
        # get groups
        try:
            grouped = self.bps.dropna().groupby(['chromosome', 'human_clust'])
        except KeyError:
            raise KeyError("Run the method `clustering_breakpoints` before"
                           " this method.")
        sorted_group = sorted(grouped.groups.items(), key=lambda x: len(x[1]),
                              reverse=True)
        cluster_df = [
            self._summarize_sub_df(cluster)
            for cluster in takewhile(lambda x: len(x[1]) >= threshold, sorted_group)
        ]
        return pd.DataFrame(
            cluster_df,
            columns=['name', 'chromosome', 'median_bpstart_human',
                     'max_end_human', 'strand_human', 'median_bpstart_virus',
                     'max_end_virus', 'number_of_read'],
        )

    def _summarize_sub_df(self, cluster):
        """ Return a summary of your cluster.
        """
        name, index = cluster
        subdf = self.bps.loc[index]
        return [
            name,
            subdf['chromosome'].max(),
            subdf['bpstart_human'].median(),
            _get_max_end(subdf, 'human'),
            subdf['strand_human'].max(),
            subdf['bpstart_virus'].median(),
            _get_max_end(subdf, 'virus'),
            len(index)
        ]

    def get_connections(self, chromosome, cluster_nb):
        """ Get connection between clusters.

        :params int cluster_nb: cluster number to find connection with other
            cluster.
        """
        bp_count = self.bps['read'].value_counts()
        # get reads name with multiple breakpoints
        multi = bp_count.loc[bp_count > 1].index
        # get reads of the cluster
        cluster_bp = self.bps.loc[
            (self.bps['chromosome'] == chromosome)
            & (self.bps['human_clust'] == cluster_nb)
        ].set_index('read')
        # get intersection reads
        intersection = cluster_bp.index.intersection(multi)
        inter_read = self.bps.loc[self.bps['read'].isin(intersection)]
        # remove the cluster of interest to get connections
        connections = inter_read.loc[inter_read['human_clust'] != cluster_nb]
        return connections

    def plot_number_bp(self, filename=None):
        """ Plot the number of reads that hold multiple breakpoint.

        :param str filename: filename to save image.

        Return image filename or show the plot in ipython.
        """
        bp_count = self.bps['read'].value_counts().value_counts()
        # Create plot
        fig, ax = init_plot("Count number of breakpoints hold by reads",
                           "Number of breakpoints", 'Count')
        ax.set_yscale('log')
        plt.bar(bp_count.index, bp_count, color="#3F5D7D")

        if filename:
            plt.savefig(filename, bbox_inches="tight",
                        facecolor=fig.get_facecolor(), edgecolor='none')
            return filename
        plt.show()


def _get_max_end(subdf, target):
    start = "bpstart_" + target
    end = "end_" + target
    if subdf[start].min() > subdf[end].min():
        return subdf[end].min()
    return subdf[end].max()


def _read_is_palindromic(bp_list):
    """ Check if a read with multiple breakpoint is palindromics.
    """
    if len(bp_list) < 2:
        return False
    bp_iter = iter(bp_list)
    prev_bp = next(bp_iter)
    for bp in bp_iter:
        if _bp_are_equal(prev_bp, bp):
            return True
        prev_bp = bp
    return False


def _bp_are_equal(bp1, bp2, margin=100):
    """ Function that determines if two breakpoints are equals.
    """
    iter_dict = zip(bp1.items(), bp2.items())
    for (k1, v1), (k2, v2) in iter_dict:
        # ignore unknown parameters
        if 'Unknown' in [v1, v2]:
            continue

        # stop parameter need special comparison
        if k1.startswith('end'):
            start = 'bpstart_' + k1.split("_")[-1]
            try:
                cond = ((bp1[start] - v1) * (bp2[start] - v2)) > 0
            except TypeError:
                cond = v1 == v2
            if cond:
                continue
            else:
                return False

        # create condition if value are string or int
        try:
            cond = abs(v1 - v2) < margin
        except TypeError:
            cond = v1 == v2
        if cond is False:
            return False
    return True
