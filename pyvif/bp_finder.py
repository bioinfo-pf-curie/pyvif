# -*- coding: utf-8 -*-

""" Breakpoint toolbox

I detect breakpoint with the minimap2 mapping and breakpoints sequences.
"""

from collections import OrderedDict

import pandas as pd
import matplotlib.pyplot as plt

from pyvif import logger, _PAF_COLNAMES
from pyvif.paftools import PAF
from pyvif.plots import init_plot, generate_plot
from pyvif.config import COLOR

_TRANSTAB = str.maketrans("ACGTacgt", "TGCAtgca")


class BreakpointFinder(object):
    """ Object that detect all breakpoint present in PacBio data.
    It uses a mapping on human and viruses genomes.
    """
    def __init__(self, human, virus):
        """.. rubric:: constructor

        :param human: PAF* format of mapping on human genome or BAM file.
        :param virus: PAF* format of mapping on virus genomes or BAM file.
        :param int threads: number of threads.

        *PAF correspond to :class:`pbcapture.bamtools.PAF` object or a PAF
        file created by :meth:`pbcapture.bamtools.PAF.to_csv`.
        """
        self.paf, self.virus_contigs = self._init_paf(human, virus)
        self.palindromics = []
        self.bps = self.find_breakpoints()
        self._clusters = None
        self.hybrid_count = self.bps['read'].unique()

    @property
    def clusters(self):
        if self._clusters is None:
            self.clustering_breakpoints()
        return self._clusters

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
        # zip function is faster than itertuples method
        iter_align = zip(
            self.paf.index, self.paf["chr"], self.paf["r_start"],
            self.paf["r_end"], self.paf["q_start"], self.paf["q_end"],
            self.paf["strand"], self.paf["mapq"], self.paf["q_length"],
        )
        # init variable as None
        prev_chr = prev_coord = prev_strand = prev_mapq = None
        read = length = None
        bp_list = list()

        for (cur_read, cur_chr, *cur_coord, cur_strand, cur_mapq,
             cur_length) in iter_align:
            # reboot analysis if read change
            if cur_read != read:
                (prev_chr, prev_coord, prev_strand, prev_mapq) = (
                    cur_chr, cur_coord, cur_strand, cur_mapq)
                # new read
                read, length = cur_read, cur_length
                continue
            elif cur_chr in self.virus_contigs:
                if prev_chr != cur_chr:
                    bp_list.append(
                        self._get_breakpoint(
                            prev_chr, prev_coord, prev_strand, prev_mapq,
                            cur_chr, cur_coord, cur_strand, read, length
                        )
                    )
            elif prev_chr in self.virus_contigs:
                bp_list.append(
                    self._get_breakpoint(
                        cur_chr, cur_coord, cur_strand, cur_mapq, prev_chr,
                        prev_coord, prev_strand, read, length
                    )
                )
            (prev_chr, prev_coord, prev_strand, prev_mapq) = (
                cur_chr, cur_coord, cur_strand, cur_mapq)
        df = pd.DataFrame(bp_list)
        logger.info("{} breakpoints are found.".format(len(bp_list)))
        return df

    def _get_breakpoint(self, human_chr, human_coord, human_strand, human_mapq,
                        virus_chr, virus_coord, virus_strand, read, length):
        """ Create breakpoint dictionnary to create dataframe.
        """
        # check if hg - hpv or hpv - hg case
        i = 0 if human_coord[2] < virus_coord[2] else 1

        # get breakpoint location
        location = virus_coord[2 + i]

        # get coord if strand is '+'
        too_far = abs(human_coord[3 - i] - virus_coord[2 + i]) > 100
        bpstart_human = end_human = "Unknown"
        if too_far:
            human_chr = "Unmapped"
            human_strand = "?"
        elif human_mapq <= 30:
            human_chr = "LowMapQ"
            human_strand = "?"
        else:
            bpstart_human = human_coord[1 - i]
            end_human = human_coord[0 + i]
        bpstart_virus = virus_coord[0 + i]
        end_virus = virus_coord[1 - i]

        # switch start / end if strand is '-'
        if human_strand == '-':
            bpstart_human, end_human = end_human, bpstart_human
        if virus_strand == '-':
            bpstart_virus, end_virus = end_virus, bpstart_virus
            # the virus strand always set as '+'
            virus_strand = '+'
            if human_strand != '?':
                human_strand = '+' if human_strand == '-' else '-'

        # build the junction plan
        virus_plan = '+>' if bpstart_virus < end_virus else '<+'
        human_plan = '{}>' if bpstart_human < end_human else '<{}'
        jct_plan = "{}/{}".format(human_plan.format(human_strand), virus_plan)

        return OrderedDict([
            ("chromosome", human_chr),
            ("bpstart_human", bpstart_human),
            ("end_human", end_human),
            ("virus_contig", virus_chr),
            ("bpstart_virus", bpstart_virus),
            ("end_virus", end_virus),
            ("jct_plan", jct_plan),
            ("read", read),
            ("length", length),
            ("location", location),
        ])

    def find_palindromics(self):
        """ Find reads that are palindromics using clustering.
        """
        logger.info("Palindromics research is running...")
        return NotImplemented

    def n50(self):
        """ Compute N50 for hybrid reads.
        """
        lengths = self.paf.groupby('q_name')\
                          .q_length.max()\
                          .sort_values(ascending=False)
        return lengths.loc[lengths.cumsum() >= sum(lengths) / 2][0]

    def clustering_breakpoints(self, human_thd=1, virus_thd=1):
        """ Generate a clustering with bpstart_human and with bpstart_virus.

        :params int human_thd: maximal distance between two human breakpoint
                               start.
        :params int virus_thd: maximal distance between two virus breakpoint
                               start.

        Add 'cluster' column in dataframe.
        """
        # clustering with human breakpoint start and chromosome
        cluster_df = self.bps[['chromosome', 'jct_plan']].copy()
        # crash if no unknown
        try:
            cond = self.bps['bpstart_human'] != 'Unknown'
            known_bp = self.bps.loc[cond]
        except TypeError:
            known_bp = self.bps
        human_cluster = known_bp.sort_values('bpstart_human')\
                                .groupby('chromosome')['bpstart_human']\
                                .diff().gt(human_thd).cumsum()
        cluster_df.loc[human_cluster.index, 'human_clust'] = human_cluster
        # create cluster with Unmapped and LowMapQ
        cluster_df['human_clust'] = cluster_df['human_clust'].fillna(
            max(human_cluster) + 1).astype(int)
        # clustering with virus breakpoint start
        virus_cluster = self.bps.sort_values('bpstart_virus')['bpstart_virus']\
                            .diff().gt(virus_thd).cumsum()
        cluster_df.loc[virus_cluster.index, 'virus_clust'] = virus_cluster
        grouped = cluster_df.groupby(['chromosome', 'human_clust',
                                      'virus_clust', 'jct_plan'])
        self._clusters = sorted(grouped.groups.items(),
                                key=lambda x: len(x[1]),
                                reverse=True)
        # add cluster rank to bps
        map_clust_to_index = {
            index: i
            for i, (__, indexes) in enumerate(self._clusters)
            for index in indexes
        }
        self.bps['cluster'] = self.bps.index.map(map_clust_to_index)

    def get_read_cluster_content(self):
        """ Return a dataframe with reads contains of clusters.
        """
        jct_count = self.bps.groupby('read').location.count()
        multiple_jct = jct_count.loc[jct_count >= 2].index
        multi_jct_df = self.bps.loc[self.bps['read'].isin(multiple_jct)]

        def group_read_info(x):
            return pd.Series({
                'count': len(x.cluster),
                'clust_content': "-".join(str(c) for c in x.cluster),
                'length': x.length.max(),
            })

        jct_content = multi_jct_df.sort_values('location')\
                                  .groupby('read')\
                                  .apply(group_read_info)
        return jct_content

    def summarise_clustering(self):
        """ Return a dataframe that summarizes bigger clusters.

        :params int threshold: minimum size of cluster.
        """
        cluster_df = [
            self._summarize_sub_df(cluster)
            for cluster in self.clusters
        ]
        return pd.DataFrame(
            cluster_df,
            columns=['rank', 'chromosome', 'median_bpstart_human',
                     'max_end_human', 'jct_plan', 'median_bpstart_virus',
                     'max_end_virus', 'number_of_read'],
        )

    def _summarize_sub_df(self, cluster):
        """ Return a summary of your cluster.
        """
        name, index = cluster
        subdf = self.bps.loc[index]
        try:
            return [
                subdf['cluster'].max(),
                subdf['chromosome'].max(),
                subdf['bpstart_human'].median(),
                _get_max_end(subdf, 'human'),
                subdf['jct_plan'].max(),
                subdf['bpstart_virus'].median(),
                _get_max_end(subdf, 'virus'),
                len(index)
            ]
        except TypeError:
            return [
                subdf['cluster'].max(),
                subdf['chromosome'].max(),
                "Unknown",
                "Unknown",
                subdf['jct_plan'].max(),
                subdf['bpstart_virus'].median(),
                _get_max_end(subdf, 'virus'),
                len(index)
            ]

    def get_bp_in_cluster(self, rank):
        """ Get breakpoint dataframe with the cluster rank.

        :params int rank: cluster rank (0 based btw).

        return dataframe with breakpoint from cluster.
        """
        cluster = self.clusters[rank]
        return self.bps.loc[cluster[1]]

    def get_bp_connections(self, rank):
        """ Get connection between clusters.

        :params int rank: cluster rank (0 based btw).

        return dataframe with connected breakpoint.
        """
        # get read name of rank
        try:
            cluster_bp = self.bps['cluster'] == rank
        except KeyError:
            self.clustering_breakpoints()
            cluster_bp = self.bps['cluster'] == rank
        self.clustering_breakpoints()
        cluster_read = self.bps.loc[cluster_bp]['read'].unique()
        connections = self.bps.loc[
            (~cluster_bp)
            & (self.bps['read'].isin(cluster_read))
        ]
        return connections

    def get_alignment_in_cluster(self, rank):
        """ Get breakpoint alignment with the cluster rank.

        :params int rank: cluster rank (0 based btw).

        return dataframe with breakpoint from cluster.
        """
        cluster = self.clusters[rank]
        reads = self.bps.loc[cluster[1]]['read'].unique()
        return self.paf.loc[reads]

    def plot_positions(self, filename=None):
        """ Barplot of reads position in reference.

        :param str filename: filename to save images.

        Return image filename.
        """
        paf = PAF(self.paf.reset_index()[_PAF_COLNAMES])
        return paf.plot_positions(filename=filename)

    def plot_bp_positions(self, filename=None):
        bp_position = self.bps['chromosome'].value_counts()
        fig, ax = init_plot(
            title="Position of breakpoints",
            xlabel="Positions",
            ylabel="Count",
            log=True
        )

        plt.xticks(fontsize=14, rotation=50)
        for tick in ax.xaxis.get_majorticklabels():
            tick.set_horizontalalignment('right')
        plt.bar(bp_position.index, bp_position, color=COLOR)
        return generate_plot(filename, fig)

    def plot_number_bp(self, filename=None):
        """ Plot the number of reads that hold multiple breakpoint.

        :param str filename: filename to save image.

        Return image filename or show the plot in ipython.
        """
        bp_count = self.bps['read'].value_counts().value_counts()
        # Create plot
        fig, ax = init_plot(
            title="Count number of breakpoints hold by reads",
            xlabel="Number of breakpoints",
            ylabel="Count",
            log=True
        )
        plt.bar(bp_count.index, bp_count, color=COLOR)
        return generate_plot(filename, fig)

    def plot_density_virus_human(self, filename=None):
        """ Plot the virus part vs the human part of reads.
        The plot shows 99% of data to zoom in.
        """
        from matplotlib import cm
        fig, ax = init_plot("", "Human Size (base)", "Virus Size (base)")

        hybrid_name = self.paf.loc[~self.paf['chr'].isin(self.virus_contigs)]\
                              .index.unique()
        hybrid = self.paf.loc[hybrid_name].reset_index()
        hybrid['size'] = hybrid['q_end'] - hybrid['q_start']
        hybrid = hybrid.set_index('q_name')

        def get_size(name):
            subdf = hybrid.loc[name]
            virus = human = 0
            for chrom, size in zip(subdf['chr'], subdf['size']):
                if chrom in self.virus_contigs:
                    virus += size
                else:
                    human += size
            return virus, human

        virus_human_list = [get_size(name) for name in hybrid_name]
        df = pd.DataFrame(virus_human_list, columns=['virus', 'human'])
        # zoom on 99% of data
        plt.hist2d(df['human'], df['virus'], bins=500, cmap=cm.viridis)
        plt.colorbar()
        ax.set_xlim(0, df['human'].quantile(0.99) + 1000)
        ax.set_ylim(0, df['virus'].quantile(0.99) + 1000)
        return generate_plot(filename, fig)

    def plot_connections_locations(self, rank, filename=None):
        """ Plot where connections are located.

        :param str filename: filename to save image.

        Return image filename or show the plot in ipython.
        """
        bp_connections = self.get_bp_connections(rank)
        posi_count = bp_connections['chromosome'].value_counts()
        fig, ax = init_plot(
            title="Position of connections",
            xlabel="Positions",
            ylabel="Count",
            log=True
        )
        plt.xticks(fontsize=14, rotation=50)
        for tick in ax.xaxis.get_majorticklabels():
            tick.set_horizontalalignment('right')
        plt.bar(posi_count.index, posi_count, color=COLOR)
        return generate_plot(filename, fig)


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
        # pass if a parameter is unknown
        if "Unknown" in [k1, v1]:
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
