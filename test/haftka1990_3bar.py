#
import pytest
import numpy as np
from main import main
#
def test():
#
    f0_ref = [5.551, 5.551, 4.441, 3.331, 2.319, 2.223, 1.933, 1.933]
    c0_ref = [1e0, 1e-6, 4e-2, 7e-2, 1e-1, 5e-9, 3e-8, 5e-7]
#
    f=main('3bar_sph')
#
    f0=[itm[0] for itm in f]
    c0=[np.maximum(np.amax(np.array(itm[1:7])),np.amax(np.absolute(np.array(itm[7:9])))) for itm in f]
#
    assert pytest.approx(f0_ref, abs=1e-3) == f0
    assert pytest.approx(c0_ref, abs=1e-3) == c0
#
if __name__ == "__main__":
    test()
#
