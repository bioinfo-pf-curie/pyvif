# -*- coding: utf-8 -*-
from pyvif import bamtools
from . import test_dir


def test_bam_to_paf():
    paf = bamtools.bam_to_paf(test_dir + '/resources/human.bam')
    assert len(paf) == 1116
    assert len(paf.columns) == 9
