"""
Pytest fixture.
"""
import sys
import os
import pytest

project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.insert(0, project_root)

@pytest.fixture(params=[1, 2, 3, 4, 5])
def integer_input(request):
    return request.param  # Returns one integer at a time from the params list
