"""
Tests functions in cypherid.py
"""

from scripts.cypherid import mock_test

def test_mock_test(integer_input):
    """
    Tests mock function.
    :param integer_input: Expected integer. Fails if cannot be converted into an int.

    """

    try:
        x = int(integer_input)
    except ValueError:
        assert False
    else:
        result = mock_test(x)
        assert result == x * 2
