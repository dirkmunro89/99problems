#
import pytest
from main import main
#
def test():
#
    f0_ref = [1.560, 1.477, 1.418, 1.383, 1.363, 1.352, 1.346, 1.343, 1.342, 1.341]
#
    f=main('5can_mma_34')
#
    f0=[itm[0] for itm in f]
    c0=[itm[1] for itm in f]
#
    assert pytest.approx(f0_ref, abs=1e-2) == f0
#
    f0_ref = [1.560,1.274,1.270,1.304,1.319,1.329,1.333,1.336,1.338,1.338,1.339,1.340,1.340]
    c0_ref = [0.000, 0.35, 0.27, 0.14, 0.08, 0.04, 0.02, 0.01, 0.008, 0.005, 0.003,0.002,0.001]
#
    f=main('5can_mma_116')
#
    f0=[itm[0] for itm in f]
    c0=[itm[1] for itm in f]
#
    assert pytest.approx(f0_ref, abs=1e-2) == f0
    assert pytest.approx(c0_ref, abs=1e-2) == c0
#
    f0_ref = [1.560,1.265,1.251,1.259,1.250,1.258,1.249,1.258,1.250,1.258,1.250,1.259,1.250,1.259]
    c0_ref = [0.000, 0.40, 0.43, 0.43, 0.44, 0.43, 0.44, 0.43,  0.44,  0.43,  0.44, 0.42, 0.44, 0.42]
#
    f=main('5can_mma_0')
#
    f0=[itm[0] for itm in f]
    c0=[itm[1] for itm in f]
#
    assert pytest.approx(f0_ref, abs=1e-2) == f0
    assert pytest.approx(c0_ref, abs=1e-2) == c0
#
    f0_ref = [1.560,1.265,1.251,1.259,1.250,1.258,1.249,1.258,1.250,1.258,1.250,1.259,1.250,1.259]
    c0_ref = [0.000, 0.40, 0.43, 0.43, 0.44, 0.43, 0.44, 0.43,  0.44,  0.43,  0.44, 0.42, 0.44, 0.42]
#
    f=main('5can_con')
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
