# -*- coding: utf-8 -*-

"""BAM toolbox"""

import pandas as pd
import pysam

from . import logger, _PAF_COLNAMES

CLIP_FLAG = {4, 5}


def bam_to_paf(filename):
    """ Read an aligned bam file and return a dataframe that corresponds to
    PAF format usable by :class:`paftools.PAF`.

    :param str filename: aligned bam file.
    """
    logger.info("Scanning input file. Please wait.")
    with pysam.AlignmentFile(filename) as bam_in:
        total_stats = [
            [
                align.reference_name,
                align.reference_start,
                align.reference_end,
                align.mapq,
                '-' if align.is_reverse else '+',
                align.query_name,
                align.infer_read_length(),
                _get_query_start(align),
                _get_query_end(align),
            ] for align in bam_in if not align.is_unmapped
        ]
    logger.info(
        "Input file contains valid {} alignments.".format(
            len(total_stats)
        )
    )
    paf = pd.DataFrame(
        total_stats,
        columns=_PAF_COLNAMES
    )
    return paf

# Attribut query_alignment_start and query_alignment_end do not resolve reverse
# mapping.
def _get_query_start(read):
    """ Get start of the query read.
    """
    strand = -1 if read.is_reverse else 0
    if read.cigar[strand][0] in CLIP_FLAG:
        return read.cigar[strand][1]
    return 0


def _get_query_end(read):
    """ Get end of the query read.
    """
    strand = 0 if read.is_reverse else -1
    if read.cigar[strand][0] in CLIP_FLAG:
        return read.infer_read_length() - read.cigar[strand][1]
    return read.infer_read_length()
