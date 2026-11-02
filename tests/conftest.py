"""Shared pytest setup for the FreePATHS test suite"""

import sys
from pathlib import Path

import pytest

# Test the working tree, not whatever freepaths happens to be installed (a non-editable
# install is a frozen snapshot of the code):
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

# freepaths.config parses sys.argv when it is first imported, and would choke on pytest's
# own arguments. Hide them, so the config falls back to default_config.py:
sys.argv = sys.argv[:1]


@pytest.fixture(autouse=True)
def restore_scattering_config():
    """Several tests point scattering.cf at a minimal stand-in config; put the real one back"""
    import freepaths.scattering as scattering
    original = scattering.cf
    yield
    scattering.cf = original
