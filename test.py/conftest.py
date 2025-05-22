import pytest
from pathlib import Path


@pytest.fixture
def data_dir():
    return Path(__file__).parent.parent / 'data'