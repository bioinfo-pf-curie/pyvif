# -*- coding: utf-8 -*-

"""Top-level package for Python Virus Integration Finder."""

__author__ = """Dimitri Desvillechabrol"""
__email__ = 'dimitri.desvillechabrol@curie.fr'
__version__ = '0.1.0'

# colnames of PAF files
_PAF_COLNAMES = ['chr', 'r_start', 'r_end', 'mapq', 'strand', 'q_name',
                 'q_length', 'q_start', 'q_end']

# add logger
import logging

handler = logging.StreamHandler()
handler.setFormatter(logging.Formatter(
    "%(asctime)s - %(levelname)-8s %(message)s",
    "%Y-%m-%d %H:%M:%S"
))
logger = logging.getLogger(__name__)
logger.setLevel('INFO')
logger.addHandler(handler)
