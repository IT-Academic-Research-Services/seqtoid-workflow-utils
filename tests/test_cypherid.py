"""
Tests functions in cypherid.py
"""

import os
from scripts.cypherid import mock_test

def test_mock_test(instr):
    """
    Tests mock function.
    :param instr: Expected integer. Fails if cannot be converted into an int.

    """

    try:
        x = int(instr)
    except ValueError:
        assert False
    else:
        result = mock_test(x)
        assert result == x * 2
