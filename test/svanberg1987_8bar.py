#
import pytest
from main import main
#
def test():
#
    f0_ref = [13.05, 12.10, 11.67, 11.65, 11.64, 11.62, 11.60, 11.56, 11.52, 11.47, \
        11.41, 11.36, 11.31, 11.24, 11.23]
#
    f=main('8bar_mma_34')
#
    f0=[itm[0] for itm in f]
#
    assert pytest.approx(f0_ref, abs=1e-2) == f0
#
    f0_ref = [13.05, 12.10, 11.67, 11.65, 11.63, 11.60, 11.53, 11.44, 11.35, 11.25, 11.23]
#
    f=main('8bar_mma_12')
#
    f0=[itm[0] for itm in f]
#
    assert pytest.approx(f0_ref, abs=1e-2) == f0
#
    f0_ref = [13.05, 12.10, 11.67, 11.65, 11.61, 11.52, 11.42, 11.28, 11.23]
#
    f=main('8bar_mma_14')
#
    f0=[itm[0] for itm in f]
#
    assert pytest.approx(f0_ref, abs=1e-2) == f0
#
#   CONLIN for this case does not correspond; it does for the 2bar problem
#
if __name__ == "__main__":
    test()
#
