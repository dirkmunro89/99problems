#
import matplotlib.pyplot as plt
#
class Prto:
#
    def __init__(self):
        self.frnt = []
        self.sgma = 1e-1
        self.gama = 1e-5
        self.beta =  1.-self.gama
#
    def pas(self,f_1,f_k,v_k,dq):
        beta=self.beta
        gama=self.gama
        sgma=self.sgma
        if f_k+gama*v_k<=min([p[0] for p in self.frnt]) or v_k<=beta*min([p[1] for p in self.frnt]):
            df = f_1 - f_k
            if df < sgma*dq and dq > 0: # approx has descended / predicted a descent, 
            #but change in actual obj is less than 10% of predicted decrease, 
            #then something not cool, we restore, and try again
                return False
            else: # either predicted an ascent (h-type), or descent in df is sufficient (f-type)
                return True
        else: return False
#
    def add(self,f_k,v_k):
        frnt=self.frnt
        frnt.append((f_k,v_k)); frnt=sorted(frnt,key=lambda x:-x[0]); tmp=[]
        for p in frnt:
            if f_k >= p[0] or v_k >= p[1]: tmp.append(p)
        self.frnt=tmp
#
    def plt(self,f_k,v_k):
        plt.clf()
        x=[itm[1] for itm in self.frnt]
        y=[itm[0] for itm in self.frnt]
        plt.grid()
        plt.plot(x, y, color='g', linewidth=0.5)
        plt.xlabel("Vio.")
        plt.ylabel("Obj.")
        plt.scatter(x,y,marker='+',s=15,linewidths=0.5,color='b')
        plt.scatter(v_k,f_k,marker='o',s=15,linewidths=0.5,color='r',zorder=10)
        plt.rcParams['axes.autolimit_mode'] = 'round_numbers'
        plt.rcParams['axes.axisbelow'] = True
        plt.savefig('prto.eps')
        plt.close()
#
