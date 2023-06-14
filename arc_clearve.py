
import matplotlib.pyplot as plt
import numpy as np
import os.path as osp
import os, sys
curpath = os.path.abspath(os.path.dirname(__file__))
sys.path.append(curpath)
import matplotlib as mpl
from matplotlib import font_manager
from sympy import *
# fontpath='/work/wangr/dxy/scripts/font/times.ttf'
fontpath='/work/wangr/dxy/scripts/font/arial.ttf'
#fontpath="D:\\JianGuoYun\\others\\keda\\old\\fermi\\font\\times.ttf"
# ------------------ font setup ----------------------#
font_properties = font_manager.FontProperties(fname=fontpath)
styles = ['normal', 'italic', 'oblique']
weights = ['light', 'normal', 'medium', 'semibold', 'bold', 'heavy', 'black']
font = {'fontproperties': font_properties,  # 'family': 'Liberation Sans'
		'style': styles[0],
		'color': 'black',
		'weight': weights[1],
		'size': 20, 
        }
plt.rc('font',family=font, size=20)
mpl.use('Agg')
############################## data processing ##############################################
class dataDriven():
    def __init__(self, filename):
        self.file_path = os.getcwd()
        self.filename = filename
        data = np.loadtxt(osp.join(self.file_path, self.filename), skiprows=0, dtype=np.float64)
        self.XX = data[:, 0]
        self.YY = data[:, 1]
        self.ZZ = data[:, 2]
        self.eq11, self.eq22, self.eq33, self.eq44, self.eq55, self.eq66 = self.get_cor()
        self.p11, self.p22, self.p33, self.p44, self.p55, self.p66 = self.intersect()
    def get_cor(self):
        x = self.XX
        y = self.YY
        xmin = np.min(x)
        xmax = np.max(x)
        ymin = np.min(y)
        ymax = np.max(y)
        b = np.array([xmax, ymin])
        c = np.array([xmin, ymax])
        yaxis = []
        for i in range(len(y)):
            if x[i] == xmin:
                yaxis.append(y[i])
        a0 = min(yaxis)
        a = np.array([xmin, a0])
        d = (b-a) + (c-a) + a
        A = (a + b)/2
        B = (c + d)/2
        C = (c + a)/2
        D = (b + d)/2
        G = (C + D)/2
        # print(np.array([A, B, C, D, G]))
        k1 = -(G[0]-C[0])/(G[1]-C[1])
        b1 = (G[1]+C[1])/2 - k1 * (G[0]+C[0])/2
        eq1 = [k1, b1]

        eq2 = [0, (G[1] + B[1])/2]

        k3 = -(d[0]-G[0])/(d[1]-G[1])
        b3 = (d[1]+G[1])/2 - k3 * (d[0]+G[0])/2
        eq3 = [k3, b3]
    
        k4 = -(G[0]-a[0])/(G[1]-a[1])
        b4 = (G[1]+a[1])/2 - k4 * (G[0]+a[0])/2
        eq4 = [k4, b4]

        eq5 = [0, (G[1]+A[1])/2]

        k6 = -(D[0]-G[0])/(D[1]-G[1])
        b6 = (D[1]+G[1])/2 - k6 * (D[0]+G[0])/2
        eq6 = [k6, b6]
        return eq1, eq2, eq3, eq4, eq5, eq6

    def intersect(self):
        eq1, eq2, eq3, eq4, eq5, eq6 = self.eq11, self.eq22, self.eq33, self.eq44, self.eq55, self.eq66
        x, y1, y2, y3, y4, y5, y6 = symbols('x y1 y2 y3 y4 y5 y6')
        y1 = eq1[0] * x + eq1[1]
        y2 = eq2[0] * x + eq2[1]
        y3 = eq3[0] * x + eq3[1]
        y4 = eq4[0] * x + eq4[1]
        y5 = eq5[0] * x + eq5[1]
        y6 = eq6[0] * x + eq6[1]
        sec1 = solve(y4-y1, x)[0]
        sec2 = solve(y2-y1, x)[0]
        sec3 = solve(y3-y2, x)[0]
        sec4 = solve(y6-y3, x)[0]
        sec5 = solve(y5-y4, x)[0]
        sec6 = solve(y6-y5, x)[0]
        p1 = [sec1, eq1[0] * sec1 + eq1[1]]
        p2 = [sec2, eq2[0] * sec2 + eq2[1]]
        p3 = [sec3, eq3[0] * sec3 + eq3[1]]
        p4 = [sec4, eq6[0] * sec4 + eq6[1]]
        p5 = [sec5, eq5[0] * sec5 + eq5[1]]
        p6 = [sec6, eq5[0] * sec6 + eq5[1]]
        return p1, p2, p3, p4, p5, p6

    def piecefunction_up(self, x1):
        p1, p2, p3, p4, p5, p6 = self.p11, self.p22, self.p33, self.p44, self.p55, self.p66
        if x1 >= p1[0] and x1 <= p2[0]:
            return self.eq11[0] * x1 + self.eq11[1]
        elif x1 > p2[0] and x1 <= p3[0]:
            return self.eq22[0] * x1 + self.eq22[1]
        elif x1 > p2[0] and x1 <= p4[0]:
            return self.eq33[0] * x1 + self.eq33[1]

    def piecefunction_dw(self, x1):
        p1, p2, p3, p4, p5, p6 = self.p11, self.p22, self.p33, self.p44, self.p55, self.p66
        if x1 >= p1[0] and x1 <= p5[0]:
            return self.eq44[0] * x1 + self.eq44[1]
        elif x1 > p5[0] and x1 <=p6[0]:
            return self.eq55[0] * x1 + self.eq55[1]
        elif x1 > p6[0] and x1 <= p4[0]:
            return self.eq66[0] * x1 + self.eq66[1]
    
    def data_pro(self):
        x_data = []
        y_data = []
        z_data = []
        for i in range(len(self.XX)):
            if self.XX[i] >= self.p11[0] and self.XX[i] <= self.p44[0]:
                if self.YY[i] >= self.piecefunction_dw(self.XX[i]) and self.YY[i] <= self.piecefunction_up(self.XX[i]):
                    x_data.append(self.XX[i])
                    y_data.append(self.YY[i])
                    z_data.append(self.ZZ[i])
        return np.array(x_data), np.array(y_data), np.array(z_data)
##################################### plot ##########################################
def plot_arc(filepath="./", filename="dos.dat_r", figname="ss_band_r.png", font_size=16):
    db = dataDriven(filename)
    x, y, z = db.data_pro()
    plt.figure(figsize=(5.0, 4.5))
    # plt.figure(figsize=(4.0, 4.5))
    ax = plt.subplot(111)
    # cm = plt.cm.jet
    # cm = plt.cm.gist_earth
    # cm = plt.cm.cubehelix
    # cm = plt.cm.seismic
    # cm = plt.cm.seismic
    # cm = plt.cm.viridis
    cm = plt.cm.hot
    # p = ax.pcolormesh(x, y, z, cmap=cm, vmin=np.min(z)+3, vmax=np.max(z)+4, shading='flat', antialiased=True)
    # p = ax.pcolormesh(x, y, z, cmap=cm, vmin=0, vmax=4, shading='flat', antialiased=True)
    # p = ax.scatter(x, y, s=5, c=z, cmap=cm, marker="o", vmin=np.min(z)+2, vmax=np.max(z)+2, alpha=0.8)
    p = ax.scatter(x, y, s=5, c=z, cmap=cm, marker="o", alpha=0.8)
    # cb = plt.colorbar(p)
    ########################################################################
    ax.set_ylim(np.min(y), np.max(y))
    ax.set_xlim(np.min(x), np.max(x))
    ax.set_xticks([])
    ax.set_yticks([])
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    plt.savefig(figname, dpi=600, bbox_inches='tight')
    plt.close()
    print("Success !")

if __name__== "__main__":
    file_path = os.getcwd()
    print("---------------------------------------------------------")
    print("begin plot ssband !")
    plot_arc(filepath=file_path, filename="arc.dat_r", figname="arc-cleave_r.png", font_size=16)
    plot_arc(filepath=file_path, filename="arc.dat_l", figname="arc-cleave_l.png", font_size=16)



