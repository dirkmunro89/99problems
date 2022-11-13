#
import pytest
from main import main
#
def test():
#
    f0_ref = [1.560, 1.336, 1.323, 1.340, 1.340, 1.340]
    c0_ref = [-.27e-16, 0.10, 0.39e-1, 0.48e-3, 0.77e-7, -0.26e-7]
#
    f=main('5can_t2r_t2d')
#
    f0=[itm[0] for itm in f]
    c0=[max(itm[1],0.) for itm in f]
#
    assert pytest.approx(f0_ref, abs=1e-3) == f0
    assert pytest.approx(c0_ref, abs=1e-1) == c0
#
    f0_ref = [1.560, 1.560, 1.331, 1.311, 1.339, 1.340, 1.340]
    c0_ref = [0e0, 2e-12, 1e-1, 7e-2, 3e-3, 6e-6, 2e-11]
#
    f=main('5can_t2r_t2e')
#
    f0=[itm[0] for itm in f]
    c0=[itm[1] for itm in f]
#
    assert pytest.approx(f0_ref, abs=1e-3) == f0
    assert pytest.approx(c0_ref, abs=1e-1) == c0
#
    f=main('5can_t2r_t2c')
#
    f0=[itm[0] for itm in f]
    c0=[max(itm[1],0.) for itm in f]
#
    assert pytest.approx(f0_ref, abs=1e-3) == f0
    assert pytest.approx(c0_ref, abs=1e-1) == c0
#
if __name__ == "__main__":
    test()
#
