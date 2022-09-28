#
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
#
#
class Enfc:
#
    def __init__(self):
        self.frnt = []
        self.sgma = 1e-1
        self.gama = 1e-4
        self.beta =  1.-self.gama
#
    def con_pas(self,g_1,g_k,q_k):
        # if feasible descent
        gama=self.gama
        if g_k[0] + gama*g_k[0]/10. < g_1[0] and max(g_k[1:])<=0.:
            return True
        # conservative
        elif all(np.greater(q_k,g_k)):
            return True
        else:
            return False
#
    def par_pas(self,f_1,f_k,v_k,p_k):
        beta=self.beta
        gama=self.gama
        sgma=self.sgma
        frnt=self.frnt
        if f_k+gama*(v_k+f_k)<min([p[0] for p in frnt]) or v_k<beta*min([p[1] for p in frnt]):
            df = f_1 - f_k # actual descent
            dq = f_1 - p_k # approximated descent
            if df < sgma*dq and dq > 0: # approx has descended / predicted a descent, 
            #but change in actual obj is less than 10% of predicted decrease, 
            #then something not cool, we restore, and try again
                return False
            else: # either predicted an ascent (h-type), or descent in df is sufficient (f-type)
                return True
        else: return False
#
    def par_add(self,f_k,v_k,k):
        frnt=self.frnt
        frnt.append((f_k,v_k,k)); frnt=sorted(frnt,key=lambda x:-x[0]); tmp=[]
        for p in frnt:
            if f_k >= p[0] or v_k >= p[1]: tmp.append(p)
        self.frnt=tmp
#
    def par_plt(self):
#
        ltx=False
        if ltx:
            matplotlib.use("pgf")
            matplotlib.rcParams.update({
                "pgf.texsystem": "pdflatex",
                'font.family': 'serif',
                'text.usetex': True,
                'pgf.rcfonts': False,
            })
#
        plt.clf()
        z=[itm[2] for itm in self.frnt]; size=list(30.*np.array(z)/(max(z)+1)+20.)
        x=[itm[1] for itm in self.frnt]
        y=[itm[0] for itm in self.frnt]
        plt.grid()
        plt.plot(x, y, color='black', linewidth=.1)
        plt.xlabel("Vio.")
        plt.ylabel("Obj.")
        plt.scatter(x,y,c=z,marker='o',s=size,linewidths=0.5,cmap='coolwarm',alpha=1.0,zorder=5)
        plt.rcParams['axes.autolimit_mode'] = 'round_numbers'
        plt.rcParams['axes.axisbelow'] = True
        plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
        plt.colorbar()
        plt.savefig('prto.eps')
        plt.close()
#
    def cnv_plt(self,h):
#
        ltx=False
        if ltx:
            matplotlib.use("pgf")
            matplotlib.rcParams.update({
                "pgf.texsystem": "pdflatex",
                'font.family': 'serif',
                'text.usetex': True,
                'pgf.rcfonts': False,
            })
            ext='.pgf'
        else:
            ext='.eps'
#
        f=[itm[0] for itm in h]
        v=[max(itm[1:]) for itm in h]
#
        plt.clf()
        plt.grid()
        plt.plot(range(len(f)), f, color='black', linewidth=1.)
        plt.xlabel("k")
        plt.ylabel("Obj.")
        plt.xlim(0,len(f))
        plt.xticks(range(len(f)))
        plt.savefig('conf'+ext)
        plt.close()
#
        plt.clf()
        plt.grid()
        plt.plot(range(len(v)), v, color='black', linewidth=1.)
        plt.xlabel("k")
        plt.ylabel("Obj.")
        plt.xlim(0,len(v))
        plt.xticks(range(len(v)))
        plt.savefig('conv'+ext)
        plt.close()
#
