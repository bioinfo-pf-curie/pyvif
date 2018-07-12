# -*- coding: utf-8 -*-

"""Console script for pyvif."""

import os
import sys

import matplotlib
matplotlib.use('Agg')
import click
from jinja2 import Environment, PackageLoader

from .bamtools import bam_to_paf
from .bp_finder import BreakpointFinder
from .paftools import PAF
from .utils import include_file, template_data, embed_png


@click.command(
    context_settings={'help_option_names': ['-h', '--help']}
)
@click.option(
    '-h', '--human',
    type=click.Path(exists=True),
    metavar='BAM',
    nargs=1,
    required=True,
    help="BAM file on human genome."
)
@click.option(
    '-v', '--virus',
    type=click.Path(exists=True),
    metavar='BAM',
    nargs=1,
    required=True,
    help="BAM file on viruses genomes."
)
@click.option(
    '-o', '--output',
    type=click.Path(),
    metavar='OUTPUT',
    nargs=1,
    default='pyvif_report.html',
    help="Report of virus integration finder."
)
@click.option(
    '-c', '--control',
    type=click.Path(exists=True),
    default=None,
    metavar='BAM',
    nargs=1,
    help="BAM file on control genes."
)
def main(human, virus, output, control):
    # do the general quality control
    virus_df = bam_to_paf(virus, add_unmapped=True)
    qc_paf = PAF(virus_df)
    # remove unmapped reads
    virus_df = virus_df.dropna()
    virus_paf = PAF(virus_df)
    # Find breakpoints
    bp_finder = BreakpointFinder(human, virus_df)
    bp_finder.clustering_breakpoints()
    summarize = bp_finder.summarize_human_clustering()

    pyvif_qc = {
        'name': os.path.basename(human).rstrip('.bam'),
        'read_count': control_paf.reads_count(),
        'length_distrib': embed_png(qc_paf.plot_length, 'filename',
                                    title="Length distribution of reads"),
        'pass_number': embed_png(qc_paf.plot_number_pass, 'filename'),
        'hpv_type': embed_png(virus_paf.plot_positions, 'filename'),
    }

    pyvif_bp = {
        'palindromic_count': len(bp_finder.palindromics),
        'cluster_table': summarize.drop('name', axis=1).to_html(
            table_id='bp_clustering',
            index=False,
            float_format='%i',
            border=0,
            classes='table table-striped table-bordered',
            justify='center'
        ),
    }

    # create the HTML report
    env = Environment(loader=PackageLoader('pyvif', 'templates'))
    # add function in global environment
    env.globals['include_file'] = include_file
    env.globals['template_data'] = template_data
    template = env.get_template('main.html')
    report_output = template.render(pyvif_qc=pyvif_qc, pyvif_bp=pyvif_bp)
    with open(output, 'w') as fout:
        print(report_output, file=fout)


def _basic_metrics_table(human_df, qc_df, control_df=None):
    total = len(qc_df['q_name'].unique())
    virus_df = qc_df.dropna()
    on_target = 
    on_target = len(virus_df['q_name'].unique())
    only_hpv = len(virus_df.loc[~virus_df['q_name'].isin(human_df['q_name'])]['q_name'].unique())
    bp_hpv = len(human_df.loc[human_df['q_name'].isin(virus_df['q_name'])]['q_name'].unique())
    print(total, on_target, only_hpv, bp_hpv)
