#
def mods(prob):
#
    if prob == '5can_t2r_spr': #tested
        from prob.p5can_t2r_spr import init
        from prob.p5can_t2r_spr import apar
        from prob.p5can_t2r_spr import simu
        from prob.p5can_t2r_spr import caml
        from subs.t2dual_spr import t2d as subs
    elif prob == '1000wm_t2r_spr':
        from prob.p1000wm_t2r_spr import init
        from prob.p1000wm_t2r_spr import apar
        from prob.p1000wm_t2r_spr import simu
        from prob.p1000wm_t2r_spr import caml
        from prob.p1000wm_t2r_spr import subs
    elif prob == '3bar_spr_t2d':
        from prob.p3bar_spr import init
        from prob.p3bar_spr import apar
        from prob.p3bar_spr import simu
        from prob.p3bar_spr import caml
        from subs.t2duel_spr import t2d as subs
    elif prob == '3bar_spr_t2c':
        from prob.p3bar_spr import init
        from prob.p3bar_spr import apar
        from prob.p3bar_spr import simu
        from prob.p3bar_spr import caml
        from subs.t2cplx_spr import t2c as subs
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
