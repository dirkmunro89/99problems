#
import pytest
from main import main
from MMA_gen import main as mma
#
def test():
#
    f_ref=mma()
    f0_ref=[itm[0] for itm in f_ref]
    g_ref=[itm[1] for itm in f_ref]
#
    f=main('Ntop_swei_enf1t')
    f0=[itm[0] for itm in f]
    g=[itm[1] for itm in f]
#
    assert pytest.approx(f0_ref, abs=1e-3) == f0
    assert pytest.approx(g_ref, abs=1e-2) == g
#
if __name__ == "__main__":
    test()
#
