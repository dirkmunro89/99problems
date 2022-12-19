#
def mods(prob):
#
    if prob == '5can_t2r_t2d':
        from prob.p5can_t2r import init
        from prob.p5can_t2r import apar
        from prob.p5can_t2r import simu
        from prob.p5can_t2r import caml
        from subs.t2dual import t2d as subs
    elif prob == '5can_t2r_t2e':
        from prob.p5can_t2r import init
        from prob.p5can_t2r import apar
        from prob.p5can_t2r import simu
        from prob.p5can_t2r import caml
        from subs.t2duel import t2d as subs
    elif prob == '5can_t2r_t2c':
        from prob.p5can_t2r import init
        from prob.p5can_t2r import apar
        from prob.p5can_t2r import simu
        from prob.p5can_t2r import caml
        from subs.t2cplx import t2c as subs
    elif prob == '1000wml_t2r_t2c':
        from prob.p1000wml_t2r import init
        from prob.p1000wml_t2r import apar
        from prob.p1000wml_t2r import simu
        from prob.p1000wml_t2r import caml
        from subs.t2cplx import t2c as subs
    elif prob == '1000wml_t2r_t2d':
        from prob.p1000wml_t2r import init
        from prob.p1000wml_t2r import apar
        from prob.p1000wml_t2r import simu
        from prob.p1000wml_t2r import caml
        from subs.t2dual import t2d as subs
    elif prob == '1000wml_t2r_t2e':
        from prob.p1000wml_t2r import init
        from prob.p1000wml_t2r import apar
        from prob.p1000wml_t2r import simu
        from prob.p1000wml_t2r import caml
        from subs.t2duel import t2d as subs
    elif prob == '3bar_sph_t2e':
        from prob.p3bar_sph import init
        from prob.p3bar_sph import apar
        from prob.p3bar_sph import simu
        from prob.p3bar_sph import caml
        from subs.t2duel import t2d as subs
    elif prob == '3bar_sph_t2c':
        from prob.p3bar_sph import init
        from prob.p3bar_sph import apar
        from prob.p3bar_sph import simu
        from prob.p3bar_sph import caml
        from subs.t2cplx import t2c as subs
    elif prob == 'Dunny_0':
        from prob.pDunny_0 import init
        from prob.pDunny_0 import apar
        from prob.pDunny_0 import simu
        from prob.pDunny_0 import caml
        from subs.t2cplx import t2c as subs
    elif prob == 'Dunny_1':
        from prob.pDunny_1 import init
        from prob.pDunny_1 import apar
        from prob.pDunny_1 import simu
        from prob.pDunny_1 import caml
        from subs.t2cplx import t2c as subs
    elif prob == 'Dunny_2':
        from prob.pDunny_2 import init
        from prob.pDunny_2 import apar
        from prob.pDunny_2 import simu
        from prob.pDunny_2 import caml
        from subs.t2milx import t2c as subs
    elif prob == 'Dunny_11':
        from prob.pDunny_11 import init
        from prob.pDunny_11 import apar
        from prob.pDunny_11 import simu
        from prob.pDunny_11 import caml
        from prob.pDunny_11 import subs
    else:
        print('Running user defined from temp.py ...')
        from temp import init
        from temp import apar
        from temp import simu
        from temp import caml
        from temp import subs
#
    return init,apar,simu,caml,subs
#
