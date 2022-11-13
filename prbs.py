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
