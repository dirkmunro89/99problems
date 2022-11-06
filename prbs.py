#
def mods(prob):
#
    if prob == '5can_mma_34': #tested
        from prob.p5can_mma_34 import init
        from prob.p5can_mma_34 import apar
        from prob.p5can_mma_34 import simu
        from prob.p5can_mma_34 import caml
        from subs.mma87dual import mma87 as subs
    elif prob == '5can_mma_116': #tested
        from prob.p5can_mma_116 import init
        from prob.p5can_mma_116 import apar
        from prob.p5can_mma_116 import simu
        from prob.p5can_mma_116 import caml
        from subs.mma87dual import mma87 as subs
    elif prob == '5can_mma_0': #tested
        from prob.p5can_mma_0 import init
        from prob.p5can_mma_0 import apar
        from prob.p5can_mma_0 import simu
        from prob.p5can_mma_0 import caml
        from subs.mma87dual import mma87 as subs
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
        from subs.mma87dual import mma87 as subs
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
        from subs.mma87dual import mma87 as subs
    elif '8bar_mma_12' in prob: #tested
        from prob.p8bar_mma_12 import init
        from prob.p8bar_mma_12 import apar
        from prob.p8bar_mma_12 import simu
        from prob.p8bar_mma_12 import caml
        from subs.mma87dual import mma87 as subs
    elif '8bar_mma_14' in prob: #tested
        from prob.p8bar_mma_14 import init
        from prob.p8bar_mma_14 import apar
        from prob.p8bar_mma_14 import simu
        from prob.p8bar_mma_14 import caml
        from subs.mma87dual import mma87 as subs
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
        from subs.mma87dual import mma87 as subs
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
#       from subs.t2dual import t2d as subs
#       from subs.t2osqp import t2osqp as subs
#       from subs.t2cvxo import t2cvxo as subs
    elif prob == 'Ntop_swei_duy8':
        from prob.pNtop_swei_duy8 import init
        from prob.pNtop_swei_duy8 import apar
        from prob.pNtop_swei_duy8 import simu
        from prob.pNtop_swei_duy8 import caml
        from prob.pNtop_swei_duy8 import subs
    elif prob == 'Ntop_swei_enf1a': # t2ar c-a 
        from prob.pNtop_swei_enf1a import init
        from prob.pNtop_swei_enf1a import apar
        from prob.pNtop_swei_enf1a import simu
        from prob.pNtop_swei_enf1a import caml
        from prob.pNtop_swei_enf1a import subs
    elif prob == 'Ntop_swei_enf1b': # t2ar t-r
        from prob.pNtop_swei_enf1b import init
        from prob.pNtop_swei_enf1b import apar
        from prob.pNtop_swei_enf1b import simu
        from prob.pNtop_swei_enf1b import caml
        from prob.pNtop_swei_enf1b import subs
    elif prob == 'Ntop_swei_enf1c': # slp aml
        from prob.pNtop_swei_enf1c import init
        from prob.pNtop_swei_enf1c import apar
        from prob.pNtop_swei_enf1c import simu
        from prob.pNtop_swei_enf1c import caml
        from prob.pNtop_swei_enf1c import subs
    elif prob == 'Ntop_swei_enf1d': # mma
        from prob.pNtop_swei_enf1d import init
        from prob.pNtop_swei_enf1d import apar
        from prob.pNtop_swei_enf1d import simu
        from prob.pNtop_swei_enf1d import caml
        from prob.pNtop_swei_enf1d import subs
    elif prob == 'Ntop_swei_enf1e': # mma
        from prob.pNtop_swei_enf1e import init
        from prob.pNtop_swei_enf1e import apar
        from prob.pNtop_swei_enf1e import simu
        from prob.pNtop_swei_enf1e import caml
        from prob.pNtop_swei_enf1e import subs
    elif prob == 'Ntop_swei_enf1t':
        from prob.pNtop_swei_enf1t import init
        from prob.pNtop_swei_enf1t import apar
        from prob.pNtop_swei_enf1t import simu
        from prob.pNtop_swei_enf1t import caml
        from prob.pNtop_swei_enf1t import subs
    elif prob == 'Ntop_swei_enf1v':
        from prob.pNtop_swei_enf1v import init
        from prob.pNtop_swei_enf1v import apar
        from prob.pNtop_swei_enf1v import simu
        from prob.pNtop_swei_enf1v import caml
        from prob.pNtop_swei_enf1v import subs
    elif prob == 'Ntop_swei_enf2a':
        from prob.pNtop_swei_enf2a import init
        from prob.pNtop_swei_enf2a import apar
        from prob.pNtop_swei_enf2a import simu
        from prob.pNtop_swei_enf2a import caml
        from prob.pNtop_swei_enf2a import subs
    elif prob == 'Ntop_swei_enf2b':
        from prob.pNtop_swei_enf2b import init
        from prob.pNtop_swei_enf2b import apar
        from prob.pNtop_swei_enf2b import simu
        from prob.pNtop_swei_enf2b import caml
        from prob.pNtop_swei_enf2b import subs
    elif prob == 'Ntop_swei_enf2c':
        from prob.pNtop_swei_enf2c import init
        from prob.pNtop_swei_enf2c import apar
        from prob.pNtop_swei_enf2c import simu
        from prob.pNtop_swei_enf2c import caml
        from prob.pNtop_swei_enf2c import subs
    elif prob == 'Ntop_swei_enf2d':
        from prob.pNtop_swei_enf2d import init
        from prob.pNtop_swei_enf2d import apar
        from prob.pNtop_swei_enf2d import simu
        from prob.pNtop_swei_enf2d import caml
        from prob.pNtop_swei_enf2d import subs
    elif prob == 'Ntop_swei_enf2e':
        from prob.pNtop_swei_enf2e import init
        from prob.pNtop_swei_enf2e import apar
        from prob.pNtop_swei_enf2e import simu
        from prob.pNtop_swei_enf2e import caml
        from prob.pNtop_swei_enf2e import subs
    elif prob == 'Ntop_mech_enf1a':
        from prob.pNtop_mech_enf1a import init
        from prob.pNtop_mech_enf1a import apar
        from prob.pNtop_mech_enf1a import simu
        from prob.pNtop_mech_enf1a import caml
        from prob.pNtop_mech_enf1a import subs
    elif prob == 'Ntop_mech_enf2a':
        from prob.pNtop_mech_enf2a import init
        from prob.pNtop_mech_enf2a import apar
        from prob.pNtop_mech_enf2a import simu
        from prob.pNtop_mech_enf2a import caml
        from prob.pNtop_mech_enf2a import subs
    elif prob == 'Ntop_mech_enf3a':
        from prob.pNtop_mech_enf3a import init
        from prob.pNtop_mech_enf3a import apar
        from prob.pNtop_mech_enf3a import simu
        from prob.pNtop_mech_enf3a import caml
        from prob.pNtop_mech_enf3a import subs
    elif prob == 'Ntop_mech_enf4a':
        from prob.pNtop_mech_enf4a import init
        from prob.pNtop_mech_enf4a import apar
        from prob.pNtop_mech_enf4a import simu
        from prob.pNtop_mech_enf4a import caml
        from prob.pNtop_mech_enf4a import subs
    elif prob == 'Ntop_sloc_enf1a':
        from prob.pNtop_sloc_enf1a import init
        from prob.pNtop_sloc_enf1a import apar
        from prob.pNtop_sloc_enf1a import simu
        from prob.pNtop_sloc_enf1a import caml
        from prob.pNtop_sloc_enf1a import subs
    elif prob == 'Thou_fle':
        from prob.pThou_fle import init
        from prob.pThou_fle import apar
        from prob.pThou_fle import simu
        from prob.pThou_fle import caml
        from prob.pThou_fle import subs
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
