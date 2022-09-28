#
def mods(prob):
#
    if prob == '5can_mma_34': #tested
        from prob.p5can_mma_34 import init
        from prob.p5can_mma_34 import apar
        from prob.p5can_mma_34 import simu
        from prob.p5can_mma_34 import caml
        from subs.mmadual import mma as subs
    elif prob == '5can_mma_116': #tested
        from prob.p5can_mma_116 import init
        from prob.p5can_mma_116 import apar
        from prob.p5can_mma_116 import simu
        from prob.p5can_mma_116 import caml
        from subs.mmadual import mma as subs
    elif prob == '5can_mma_0': #tested
        from prob.p5can_mma_0 import init
        from prob.p5can_mma_0 import apar
        from prob.p5can_mma_0 import simu
        from prob.p5can_mma_0 import caml
        from subs.mmadual import mma as subs
    elif prob == '5can_t2r': #tested
        from prob.p5can_t2r import init
        from prob.p5can_t2r import apar
        from prob.p5can_t2r import simu
        from prob.p5can_t2r import caml
        from subs.t2dual import t2d as subs
    elif prob == '5can_con': #tested
        from prob.p5can_mma_0 import init
        from prob.p5can_mma_0 import apar
        from prob.p5can_mma_0 import simu
        from prob.p5can_mma_0 import caml
        from subs.condual import con as subs
    elif prob == '2bar_mma': #tested
        from prob.p2bar_mma import init
        from prob.p2bar_mma import apar
        from prob.p2bar_mma import simu
        from prob.p2bar_mma import caml
        from subs.mmadual import mma as subs
    elif prob == '2bar_con': #tested
        from prob.p2bar import init
        from prob.p2bar import apar
        from prob.p2bar import simu
        from prob.p2bar import caml
        from subs.condual import con as subs
    elif prob == '2bar_slp': #tested
        from prob.p2bar import init
        from prob.p2bar import apar
        from prob.p2bar import simu
        from prob.p2bar import caml
        from subs.t2dual import t2d as subs
    elif '8bar_mma_34' in prob: #tested
        from prob.p8bar_mma_34 import init
        from prob.p8bar_mma_34 import apar
        from prob.p8bar_mma_34 import simu
        from prob.p8bar_mma_34 import caml
        from subs.mmadual import mma as subs
    elif '8bar_mma_12' in prob: #tested
        from prob.p8bar_mma_12 import init
        from prob.p8bar_mma_12 import apar
        from prob.p8bar_mma_12 import simu
        from prob.p8bar_mma_12 import caml
        from subs.mmadual import mma as subs
    elif '8bar_mma_14' in prob: #tested
        from prob.p8bar_mma_14 import init
        from prob.p8bar_mma_14 import apar
        from prob.p8bar_mma_14 import simu
        from prob.p8bar_mma_14 import caml
        from subs.mmadual import mma as subs
    elif prob == '10bar_t2r':
        from prob.p10bar_t2r import init
        from prob.p10bar_t2r import apar
        from prob.p10bar_t2r import simu
        from prob.p10bar_t2r import caml
        from subs.t2dual import t2d as subs
    elif prob == '10bar_t2c':
        from prob.p10bar_t2c import init
        from prob.p10bar_t2c import apar
        from prob.p10bar_t2c import simu
        from prob.p10bar_t2c import caml
        from subs.t2dual import t2d as subs
    elif '8bar_con' in prob:
        from prob.p8bar import init
        from prob.p8bar import apar
        from prob.p8bar import simu
        from prob.p8bar import caml
        from subs.condual import con as subs
    elif prob == '5can_mma': #tested
        from prob.p5can_mma import init
        from prob.p5can_mma import apar
        from prob.p5can_mma import simu
        from prob.p5can_mma import caml
        from subs.mmadual import mma as subs
    elif prob == '5can_t2m':
        from prob.p5can_t2m import init
        from prob.p5can_t2m import apar
        from prob.p5can_t2m import simu
        from prob.p5can_t2m import caml
        from subs.t2dual import t2d as subs
    elif prob == 'arch':
        from arch import init
        from arch import simu
        from arch import caml
#       from qlpml_arch import qlp as subs
#       from dlpml_arch import dlp as subs
    elif prob == 'beam':
        from beam import init
        from beam import simu
        from beam import caml
    elif prob == 'Ntop_eoc':
        from prob.pNtop_eoc import init
        from prob.pNtop_eoc import apar
        from prob.pNtop_eoc import simu
        from prob.pNtop_eoc import caml
        from subs.eocdual import eoc as subs
    elif prob == 'Ntop_swei_duy8':
        from prob.pNtop_swei_duy8 import init
        from prob.pNtop_swei_duy8 import apar
        from prob.pNtop_swei_duy8 import simu
        from prob.pNtop_swei_duy8 import caml
        from prob.pNtop_swei_duy8 import subs
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
def keep(prob):
#
    if prob == 'arch': 
        m_r=1e8; m_a=0.2
        a_f=[0.7,1.2]
        kmax=2000; conv=1e-2
    elif prob == 'comp': 
        m_r=1e8; m_a=0.2
        a_f=[0.7,1.2]
        kmax=2000; conv=1e-2
    elif prob == 'beam': 
        m_r=1e8; m_a=0.2
        a_f=[0.7,1.1]
        kmax=1000; conv=1e-2
#
    return m_r, m_a, a_f, kmax, conv
#
