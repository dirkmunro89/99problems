#
import pytest
from main import main
from topopt_cholmod import main as topchl
#
def test():
#
    nelx=180
    nely=60
    volfrac=0.4
    rmin=5.4
    penal=3.0
    ft=1 # ft==0 -> sens, ft==1 -> dens
#
    f_ref=topchl(nelx,nely,volfrac,penal,rmin,ft)
    f0_ref=[itm[0] for itm in f_ref]
    g_ref=[itm[1] for itm in f_ref]
#
    f=main('Ntop_eoc')
    f0=[itm[0] for itm in f]
    g=[itm[1] for itm in f]
#
    assert pytest.approx(f0_ref, abs=1e-1) == f0
    assert pytest.approx(g_ref, abs=1e-3) == g
#
    f=main('Ntop_eoc_el')
    f0=[itm[0] for itm in f]
    g=[itm[1] for itm in f]
#
    assert pytest.approx(f0_ref, abs=1e-1) == f0
    assert pytest.approx(g_ref, abs=1e-3) == g
#
if __name__ == "__main__":
    test()
#
