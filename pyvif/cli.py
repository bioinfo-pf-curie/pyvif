# -*- coding: utf-8 -*-

"""Console script for pyvif."""

import os
import sys

import click
from jinja2 import Environment, PackageLoader

from .bamtools import bam_to_paf
from .bp_finder import BreakpointFinder
from .paftools import PAF
from .utils import include_file, template_data, embed_png


#@click.command(
#    context_settings={'help_option_names': ['-h', '--help']}
#)
#@click.option(
#    '-h', '--human',
#    type=click.Path(exists=True),
#    metavar='BAM',
#    nargs=1,
#    required=True,
#    help="BAM file on human genome."
#)
#@click.option(
#    '-v', '--virus',
#    type=click.Path(exists=True),
#    metavar='BAM',
#    nargs=1,
#    required=True,
#    help="BAM file on viruses genomes."
#)
#@click.option(
#    '-o', '--output',
#    type=click.Path(),
#    metavar='OUTPUT',
#    nargs=1,
#    default='pyvif_report.html',
#    help="Report of virus integration finder."
#)
def main(human, virus, output):
    # do the general quality control
    virus = bam_to_paf(virus, add_unmapped=True)
    control_paf = PAF(virus)
    pyvif_qc = {
        "name": os.path.basename(human).rstrip('.bam'),
        "read_count": control_paf.reads_count(),
        "length_distrib": embed_png(control_paf.plot_length, 'filename',
                                    title="Length distribution of reads"),
    }

    # Find breakpoints
    bp_finder = BreakpointFinder(human, virus)

    # create the HTML report
    env = Environment(loader=PackageLoader('pyvif', 'templates'))
    # add function in global environment
    env.globals['include_file'] = include_file
    env.globals['template_data'] = template_data
    template = env.get_template('main.html')
    report_output = template.render(data_info=data_info)
    with open(output, 'w') as fout:
        print(report_output, file=fout)
