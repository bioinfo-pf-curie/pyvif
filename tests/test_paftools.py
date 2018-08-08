# -*- coding: utf-8 -*-
import pytest

from pyvif import paftools
from . import test_dir


@pytest.fixture
def test_paf():
    # init with bamfile
    paf1 = paftools.PAF(test_dir + "/resources/human.bam")
    # init with df
    paf2 = paftools.PAF(paf1.df)

    assert len(paf1.df) == 1116
    assert paf1.reads_count() == paf2.reads_count()
    return paf1

def test_paf_plot_length(test_paf, tmpdir):
    test_paf.plot_length(
        filename=tmpdir.join("plot_length.png")
    )

def test_paf_plot_positions(test_paf, tmpdir):
    test_paf.plot_positions(
        filename=tmpdir.join("plot_position.png")
    )

def test_paf_plot_length_vs_mapping(test_paf, tmpdir):
    test_paf.plot_length_vs_mapping(
        filename=tmpdir.join("plot_length_vs_mapping.png")
    )

def test_paf_plot_number_of_pass(test_paf, tmpdir):
    test_paf.plot_number_pass(
        filename=tmpdir.join("plot_number_pass.png")
    )
