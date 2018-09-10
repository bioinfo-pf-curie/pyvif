# -*- coding: utf-8 -*-
import pytest

from pyvif import bp_finder
from . import test_dir


@pytest.fixture
def test_bp_finder():
    bp_found = bp_finder.BreakpointFinder(
        human=test_dir + "/resources/human.bam",
        virus=test_dir + "/resources/virus.bam"
    )
    assert len(bp_found.bps) == 1079
    assert len(bp_found.bps.columns) == 9
    return bp_found

def test_bp_clustering(test_bp_finder):
    test_bp_finder.clustering_breakpoints()
    assert len(test_bp_finder.bps.columns) == 10
    summarise = test_bp_finder.summarise_clustering()
    assert len(summarise) == 13
    assert summarise.loc[0,"number_of_read"] == 1042
    assert len(test_bp_finder.get_bp_in_cluster(0)) == 1042
    assert len(test_bp_finder.get_bp_connections(0)) == 15
    assert len(test_bp_finder.get_alignment_in_cluster(0)) == 2323

def test_plot_number_bp(test_bp_finder, tmpdir):
    test_bp_finder.plot_number_bp(tmpdir.join("plot_number_bp.png"))

def test_plot_connections_locations(test_bp_finder, tmpdir):
    test_bp_finder.plot_connections_locations(0, tmpdir.join("plot_connections.png"))

def test_plot_density_virus_human(test_bp_finder, tmpdir):
    test_bp_finder.plot_density_virus_human(tmpdir.join("plot_density.png"))
