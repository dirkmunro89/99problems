#
import pytest
import numpy as np
from main import main
#
def test():
#
    f0_ref=950.
#
    f=main('1000wml_t2r_t2d')
#
    assert pytest.approx(f0_ref, abs=1e-3) == f[-1][0]
    assert pytest.approx(36, abs=0.) == len(f)
#
    f=main('1000wml_t2r_t2e')
#
    assert pytest.approx(f0_ref, abs=1e-3) == f[-1][0]
    assert pytest.approx(40, abs=0.) == len(f)
#
    f=main('1000wml_t2r_t2c')
#
    assert pytest.approx(f0_ref, abs=1e-3) == f[-1][0]
    assert pytest.approx(24, abs=0.) == len(f)
#
if __name__ == "__main__":
    test()
#
