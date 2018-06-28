# -*- coding: utf-8 -*-

"""Console script for pyvif."""

import sys
import click
from jinja2 import Environment, PackageLoader

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
    '-v', '--viruses',
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
def main(human, viruses, output):
   pass 
