#
import pytest
from main import main
#
def test():
#
    f0_ref = [1.560, 1.336, 1.323, 1.340, 1.340, 1.340]
    c0_ref = [-.27e-16, 0.10, 0.39e-1, 0.48e-3, 0.77e-7, -0.26e-7]
#
    f=main('5can_t2r')
#
    f0=[itm[0] for itm in f]
    c0=[itm[1] for itm in f]
#
    assert pytest.approx(f0_ref, abs=1e-3) == f0
    assert pytest.approx(c0_ref, abs=1e-2) == c0
#
    f0_ref = [1.560,1.446,1.373,1.342,1.339,1.339,1.339,1.339,1.339,1.339,1.339,1.339,1.339,1.339]
    c0_ref = [-.27e-6,-.18e-1,-0.69e-2,0.44e-2,0.90e-3,0.15e-3,0.41e-4,0.75e-5,0.77e-6,0.53e-7,\
        0.28e-7,0.71e-8,0.64e-7,0.18e-9]
#
    f=main('5can_mma')
#
    f0=[itm[0] for itm in f]
    c0=[itm[1] for itm in f]
#
    assert pytest.approx(f0_ref, abs=1e-2) == f0
    assert pytest.approx(c0_ref, abs=1e-2) == c0
#
if __name__ == "__main__":
    test()
#
