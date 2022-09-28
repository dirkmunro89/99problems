#
import pytest
from main import main
#
def test():
#
    f0_ref = [1.68, 1.42, 1.27, 1.38, 1.28, 1.38, 1.28, 1.38]
    s1_ref = [0.92, 1.11, 1.22, 1.14, 1.21, 1.14, 1.21, 1.14]
#
    f=main('2bar_slp')
#
    f0=[itm[0] for itm in f]
    s1=[itm[1]+1 for itm in f]
#
    assert pytest.approx(f0_ref, abs=1e-2) == f0
    assert pytest.approx(s1_ref, abs=1e-2) == s1
#
    f0_ref = [1.68, 1.43, 1.37, 1.44, 1.47, 1.51]
    s1_ref = [0.92, 1.10, 1.13, 1.10, 1.03, 1.00]
#
    f=main('2bar_mma')
#
    f0=[itm[0] for itm in f]
    s1=[itm[1]+1 for itm in f]
#
    assert pytest.approx(f0_ref, abs=1e-2) == f0
    assert pytest.approx(s1_ref, abs=1e-2) == s1
#
    f0_ref = [1.68, 1.43, 1.49, 1.43, 1.49, 1.43, 1.49, 1.43]
    s1_ref = [0.92, 1.11, 1.04, 1.11, 1.04, 1.11, 1.04, 1.11]
#
    f=main('2bar_con')
#
    f0=[itm[0] for itm in f]
    s1=[itm[1]+1 for itm in f]
#
    assert pytest.approx(f0_ref, abs=1e-2) == f0
    assert pytest.approx(s1_ref, abs=1e-2) == s1
#
if __name__ == "__main__":
    test()
#
