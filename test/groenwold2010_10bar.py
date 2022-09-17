#
import pytest
from main import main
#
def test():
#
    f0_ref=[2098,2388.6,1630.4,1766.4,1587.2,1511.4,1506.6,1506.0,1505.8,1505.5,1505.2,1505.0,1504.7]
#
    f=main('10bar_t2r')
#
    f0=[itm[0] for itm in f]
#
    assert pytest.approx(f0_ref, abs=1e0) == f0
#
    f0_ref = [2098,2310.1,1510.4,1705.4,1510.0,1497.5,1497.5,1497.5]
#
    f=main('10bar_t2c')
#
    f0=[itm[0] for itm in f]
#
    assert pytest.approx(f0_ref, abs=2e0) == f0
#
if __name__ == "__main__":
    test()
#
