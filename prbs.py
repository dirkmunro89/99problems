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
    elif prob == 'Ntop_eoc':
        from prob.pNtop_eoc import init
        from prob.pNtop_eoc import apar
        from prob.pNtop_eoc import simu
        from prob.pNtop_eoc import caml
        from subs.eocdual import eoc as subs
    elif prob == 'Ntop_eoc_el':
        from prob.pNtop_eoc_el import init
        from prob.pNtop_eoc_el import apar
        from prob.pNtop_eoc_el import simu
        from prob.pNtop_eoc_el import caml
        from subs.eocdual import eoc as subs
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
    elif prob == 'Ntop_swei_enf1aa': # t2ar c-a 
        from prob.pNtop_swei_enf1aa import init
        from prob.pNtop_swei_enf1aa import apar
        from prob.pNtop_swei_enf1aa import simu
        from prob.pNtop_swei_enf1aa import caml
        from prob.pNtop_swei_enf1aa import subs
    elif prob == 'Ntop_swei_enf1b': # t2ar t-r
        from prob.pNtop_swei_enf1b import init
        from prob.pNtop_swei_enf1b import apar
        from prob.pNtop_swei_enf1b import simu
        from prob.pNtop_swei_enf1b import caml
        from prob.pNtop_swei_enf1b import subs
    elif prob == 'Ntop_swei_enf1bb': # t2ar t-r
        from prob.pNtop_swei_enf1bb import init
        from prob.pNtop_swei_enf1bb import apar
        from prob.pNtop_swei_enf1bb import simu
        from prob.pNtop_swei_enf1bb import caml
        from prob.pNtop_swei_enf1bb import subs
    elif prob == 'Ntop_swei_enf1c': # slp aml
        from prob.pNtop_swei_enf1c import init
        from prob.pNtop_swei_enf1c import apar
        from prob.pNtop_swei_enf1c import simu
        from prob.pNtop_swei_enf1c import caml
        from prob.pNtop_swei_enf1c import subs
    elif prob == 'Ntop_swei_enf1cc': # slp aml
        from prob.pNtop_swei_enf1cc import init
        from prob.pNtop_swei_enf1cc import apar
        from prob.pNtop_swei_enf1cc import simu
        from prob.pNtop_swei_enf1cc import caml
        from prob.pNtop_swei_enf1cc import subs
    elif prob == 'Ntop_swei_enf1d': # mma
        from prob.pNtop_swei_enf1d import init
        from prob.pNtop_swei_enf1d import apar
        from prob.pNtop_swei_enf1d import simu
        from prob.pNtop_swei_enf1d import caml
        from prob.pNtop_swei_enf1d import subs
    elif prob == 'Ntop_swei_enf1dd': # mma
        from prob.pNtop_swei_enf1dd import init
        from prob.pNtop_swei_enf1dd import apar
        from prob.pNtop_swei_enf1dd import simu
        from prob.pNtop_swei_enf1dd import caml
        from prob.pNtop_swei_enf1dd import subs
    elif prob == 'Ntop_swei_enf1e': # gmma
        from prob.pNtop_swei_enf1e import init
        from prob.pNtop_swei_enf1e import apar
        from prob.pNtop_swei_enf1e import simu
        from prob.pNtop_swei_enf1e import caml
        from prob.pNtop_swei_enf1e import subs
    elif prob == 'Ntop_swei_enf1ee': # gmma
        from prob.pNtop_swei_enf1ee import init
        from prob.pNtop_swei_enf1ee import apar
        from prob.pNtop_swei_enf1ee import simu
        from prob.pNtop_swei_enf1ee import caml
        from prob.pNtop_swei_enf1ee import subs
    elif prob == 'Ntop_swei_enf1f': # lpt2r
        from prob.pNtop_swei_enf1f import init
        from prob.pNtop_swei_enf1f import apar
        from prob.pNtop_swei_enf1f import simu
        from prob.pNtop_swei_enf1f import caml
        from prob.pNtop_swei_enf1f import subs
    elif prob == 'Ntop_swei_enf1ff': # lpt2r
        from prob.pNtop_swei_enf1ff import init
        from prob.pNtop_swei_enf1ff import apar
        from prob.pNtop_swei_enf1ff import simu
        from prob.pNtop_swei_enf1ff import caml
        from prob.pNtop_swei_enf1ff import subs
    elif prob == 'Ntop_swei_enf1g': # lpt2r
        from prob.pNtop_swei_enf1g import init
        from prob.pNtop_swei_enf1g import apar
        from prob.pNtop_swei_enf1g import simu
        from prob.pNtop_swei_enf1g import caml
        from prob.pNtop_swei_enf1g import subs
    elif prob == 'Ntop_swei_enf1gg': # lpt2r
        from prob.pNtop_swei_enf1gg import init
        from prob.pNtop_swei_enf1gg import apar
        from prob.pNtop_swei_enf1gg import simu
        from prob.pNtop_swei_enf1gg import caml
        from prob.pNtop_swei_enf1gg import subs
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
    elif prob == 'Ntop_swei_enf1t_el':
        from prob.pNtop_swei_enf1t_el import init
        from prob.pNtop_swei_enf1t_el import apar
        from prob.pNtop_swei_enf1t_el import simu
        from prob.pNtop_swei_enf1t_el import caml
        from prob.pNtop_swei_enf1t_el import subs
    elif prob == 'Ntop_swei_enf1v_el':
        from prob.pNtop_swei_enf1v_el import init
        from prob.pNtop_swei_enf1v_el import apar
        from prob.pNtop_swei_enf1v_el import simu
        from prob.pNtop_swei_enf1v_el import caml
        from prob.pNtop_swei_enf1v_el import subs
    elif prob == 'Ntop_swei_enf2a':
        from prob.pNtop_swei_enf2a import init
        from prob.pNtop_swei_enf2a import apar
        from prob.pNtop_swei_enf2a import simu
        from prob.pNtop_swei_enf2a import caml
        from prob.pNtop_swei_enf2a import subs
    elif prob == 'Ntop_swei_enf2aa':
        from prob.pNtop_swei_enf2aa import init
        from prob.pNtop_swei_enf2aa import apar
        from prob.pNtop_swei_enf2aa import simu
        from prob.pNtop_swei_enf2aa import caml
        from prob.pNtop_swei_enf2aa import subs
    elif prob == 'Ntop_swei_enf2b':
        from prob.pNtop_swei_enf2b import init
        from prob.pNtop_swei_enf2b import apar
        from prob.pNtop_swei_enf2b import simu
        from prob.pNtop_swei_enf2b import caml
        from prob.pNtop_swei_enf2b import subs
    elif prob == 'Ntop_swei_enf2bb':
        from prob.pNtop_swei_enf2bb import init
        from prob.pNtop_swei_enf2bb import apar
        from prob.pNtop_swei_enf2bb import simu
        from prob.pNtop_swei_enf2bb import caml
        from prob.pNtop_swei_enf2bb import subs
    elif prob == 'Ntop_swei_enf2c':
        from prob.pNtop_swei_enf2c import init
        from prob.pNtop_swei_enf2c import apar
        from prob.pNtop_swei_enf2c import simu
        from prob.pNtop_swei_enf2c import caml
        from prob.pNtop_swei_enf2c import subs
    elif prob == 'Ntop_swei_enf2cc':
        from prob.pNtop_swei_enf2cc import init
        from prob.pNtop_swei_enf2cc import apar
        from prob.pNtop_swei_enf2cc import simu
        from prob.pNtop_swei_enf2cc import caml
        from prob.pNtop_swei_enf2cc import subs
    elif prob == 'Ntop_swei_enf2d':
        from prob.pNtop_swei_enf2d import init
        from prob.pNtop_swei_enf2d import apar
        from prob.pNtop_swei_enf2d import simu
        from prob.pNtop_swei_enf2d import caml
        from prob.pNtop_swei_enf2d import subs
    elif prob == 'Ntop_swei_enf2dd':
        from prob.pNtop_swei_enf2dd import init
        from prob.pNtop_swei_enf2dd import apar
        from prob.pNtop_swei_enf2dd import simu
        from prob.pNtop_swei_enf2dd import caml
        from prob.pNtop_swei_enf2dd import subs
    elif prob == 'Ntop_swei_enf2e':
        from prob.pNtop_swei_enf2e import init
        from prob.pNtop_swei_enf2e import apar
        from prob.pNtop_swei_enf2e import simu
        from prob.pNtop_swei_enf2e import caml
        from prob.pNtop_swei_enf2e import subs
    elif prob == 'Ntop_swei_enf2ee':
        from prob.pNtop_swei_enf2ee import init
        from prob.pNtop_swei_enf2ee import apar
        from prob.pNtop_swei_enf2ee import simu
        from prob.pNtop_swei_enf2ee import caml
        from prob.pNtop_swei_enf2ee import subs
    elif prob == 'Ntop_swei_enf2f':
        from prob.pNtop_swei_enf2f import init
        from prob.pNtop_swei_enf2f import apar
        from prob.pNtop_swei_enf2f import simu
        from prob.pNtop_swei_enf2f import caml
        from prob.pNtop_swei_enf2f import subs
    elif prob == 'Ntop_swei_enf2ff':
        from prob.pNtop_swei_enf2ff import init
        from prob.pNtop_swei_enf2ff import apar
        from prob.pNtop_swei_enf2ff import simu
        from prob.pNtop_swei_enf2ff import caml
        from prob.pNtop_swei_enf2ff import subs
    elif prob == 'Ntop_swei_enf2g':
        from prob.pNtop_swei_enf2g import init
        from prob.pNtop_swei_enf2g import apar
        from prob.pNtop_swei_enf2g import simu
        from prob.pNtop_swei_enf2g import caml
        from prob.pNtop_swei_enf2g import subs
    elif prob == 'Ntop_swei_enf2gg':
        from prob.pNtop_swei_enf2gg import init
        from prob.pNtop_swei_enf2gg import apar
        from prob.pNtop_swei_enf2gg import simu
        from prob.pNtop_swei_enf2gg import caml
        from prob.pNtop_swei_enf2gg import subs
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
    elif prob == 'Ntop_thme_enf1a': # 
        from prob.pNtop_thme_enf1a import init
        from prob.pNtop_thme_enf1a import apar
        from prob.pNtop_thme_enf1a import simu
        from prob.pNtop_thme_enf1a import caml
        from prob.pNtop_thme_enf1a import subs
    elif prob == 'Ntop_thme_enf1b': # 
        from prob.pNtop_thme_enf1b import init
        from prob.pNtop_thme_enf1b import apar
        from prob.pNtop_thme_enf1b import simu
        from prob.pNtop_thme_enf1b import caml
        from prob.pNtop_thme_enf1b import subs
    elif prob == 'Ntop_thme_enf1c': # 
        from prob.pNtop_thme_enf1c import init
        from prob.pNtop_thme_enf1c import apar
        from prob.pNtop_thme_enf1c import simu
        from prob.pNtop_thme_enf1c import caml
        from prob.pNtop_thme_enf1c import subs
    elif prob == 'Ntop_thme_enf1d': # 
        from prob.pNtop_thme_enf1d import init
        from prob.pNtop_thme_enf1d import apar
        from prob.pNtop_thme_enf1d import simu
        from prob.pNtop_thme_enf1d import caml
        from prob.pNtop_thme_enf1d import subs
    elif prob == 'Ntop_thme_enf1e': # 
        from prob.pNtop_thme_enf1e import init
        from prob.pNtop_thme_enf1e import apar
        from prob.pNtop_thme_enf1e import simu
        from prob.pNtop_thme_enf1e import caml
        from prob.pNtop_thme_enf1e import subs
    elif prob == 'Ntop_thme_enf2a': # t2ar c-a 
        from prob.pNtop_thme_enf2a import init
        from prob.pNtop_thme_enf2a import apar
        from prob.pNtop_thme_enf2a import simu
        from prob.pNtop_thme_enf2a import caml
        from prob.pNtop_thme_enf2a import subs
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
    elif prob == '3bar_sph':
        from prob.p3bar_sph import init
        from prob.p3bar_sph import apar
        from prob.p3bar_sph import simu
        from prob.p3bar_sph import caml
        from prob.p3bar_sph import subs
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
