# -*- coding: utf-8 -*-

"""Console script for pyvif."""

import os

import click
from jinja2 import Environment, PackageLoader

# change backend of matplotlib
import matplotlib
matplotlib.use('Agg')
# change plot background color
from pyvif import config
config.BACKGROUND = "#fafafa"

from pyvif.bamtools import bam_to_paf
from pyvif.bp_finder import BreakpointFinder
from pyvif.paftools import PAF
from pyvif.utils import include_file, template_data, embed_png


@click.command(
    context_settings={'help_option_names': ['-h', '--help']}
)
@click.option(
    '-u', '--human',
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
    qc_df = bam_to_paf(virus, add_unmapped=True)
    qc_paf = PAF(qc_df)
    # get reads that hit control genes
    control_df = bam_to_paf(control) if control is not None else None
    # remove unmapped reads
    virus_df = qc_df.dropna()
    virus_paf = PAF(virus_df)
    # Find breakpoints
    human_df = bam_to_paf(human)
    bp_finder = BreakpointFinder(human_df, virus_df)
    bp_finder.clustering_breakpoints()
    summarize = bp_finder.summarize_human_clustering()

    pyvif_qc = {
        'name': os.path.basename(human).rstrip('.bam'),
        'basic_metrics': _basic_metrics_table(human_df, qc_df, control_df),
        'length_distrib': embed_png(qc_paf.plot_length, 'filename',
                                    title="Read length distribution"),
        'pass_number': embed_png(qc_paf.plot_number_pass, 'filename'),
        'hpv_type': embed_png(virus_paf.plot_positions, 'filename'),
    }

    pyvif_bp = {
        'palindromic_count': len(bp_finder.palindromics),
        'bp_count': embed_png(bp_finder.plot_number_bp, 'filename'),
        'cluster_table': summarize.drop('name', axis=1).to_html(
            table_id='bp_clustering',
            index=False,
            float_format=lambda x: "{:.0f}".format(x),
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
    # count subreads
    virus_df = qc_df.dropna()
    total = len(qc_df['q_name'].unique())
    basic_metrics = {
        'Total': total,
        "On target": len(virus_df['q_name'].unique()), 
        "On virus": len(virus_df['q_name'].unique()),
        "Only virus": len(
            virus_df.loc[~virus_df['q_name'].isin(human_df['q_name'])]['q_name'].unique()
        ),
        "Breakpoint virus": len(
            human_df.loc[human_df['q_name'].isin(virus_df['q_name'])]['q_name'].unique()
        ),
    }
    # add control count
    try:
        on_control = len(control_df['q_name'].unique())
        basic_metrics["On control"] = on_control
        basic_metrics["On target"] += on_control
    except TypeError:
        pass
    # count percentage
    basic_metrics = {
        k: (v, "{:.2f}%".format(v / total * 100)) for k, v in basic_metrics.items()
    }
    return basic_metrics
