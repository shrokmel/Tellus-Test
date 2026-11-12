import os, shutil, sys
import pytest

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, os.path.join(ROOT, "sim"))
sys.path.insert(0, os.path.join(ROOT, "organizer"))

needs_cc = pytest.mark.skipif(shutil.which("cc") is None, reason="needs a C compiler (cc) to build sim/sky.c")


@pytest.fixture(scope="session")
def root():
    return ROOT
