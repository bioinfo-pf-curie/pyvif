# -*- coding: utf-8 -*-
import click.testing as cli_test
import pytest
import os

from pyvif import cli
from . import test_dir


@pytest.fixture(scope="module")
def runner():
    return cli_test.CliRunner()

def test_pyvif(runner, tmpdir):
    result = runner.invoke(cli.main,[
        "--human", test_dir + "/resources/human.bam",
        "--virus", test_dir + "/resources/virus.bam",
        "--output", tmpdir.join("test_pyvif.html"),
    ])
    assert result.exit_code == 0
    assert os.path.exists(tmpdir.join("test_pyvif.html"))

def test_default_pyvif(runner):
    result = runner.invoke(cli.main, input='\n')
    assert result.exit_code == 2
