"""
test_main.py -- smoke test
"""

import os
import pytest

from click.testing import CliRunner

from feps import main


@pytest.fixture(scope='module')
def runner():
    return CliRunner()

def test_main(runner):
    result = runner.invoke(main.cli, ['-p', 'RT', '-m', 'LY12'])
    assert result.exit_code == 0
