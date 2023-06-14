# -*- coding: utf-8 -*-
"""
Created on 19:38 07-11 2022

@author: XY Ding
mail to: dxy_vasp@163.com
python3: arc.py
"""
import numpy as np
from matplotlib import pyplot as plt
import os, sys
import os.path as osp
import matplotlib as mpl
mpl.use("agg")
curpath = os.path.abspath(os.path.dirname(__file__))
sys.path.append(curpath)
import func.funcs as F

##################### process wannier tools ###########################
class ribbon():
    def __init__(self):
        self.soft_path = curpath
        self.file_path = os.getcwd()
        self.Nk1 = int(2*F.readNk1_WtOUT(self.file_path, "WT.out"))

    def data_process(self, colorZ1):
        colorZ = (colorZ1 - np.average(colorZ1))*-1
        data = colorZ/np.max(colorZ)
        color_data = 1/(1+np.exp(-5*data))
        print(np.min(color_data))
        print(np.max(color_data))
        print(np.average(color_data))
        for ii in range(len(color_data)):
            if 0.4 < color_data[ii] <= 0.5:
                color_data[ii] = 0.4
            elif 0.5 < color_data[ii] <= 0.7:
                color_data[ii] = 0.7
        return color_data

    def plot_ribbon(self, figs, eng):
        data = np.loadtxt(osp.join(self.file_path, "ribbonek.dat"), skiprows=0, dtype=float)
        rgbv = self.data_process(data[:, 2])
        plt.figure(figsize=figs)
        ax = plt.subplot(111)
        print("begin plot!")
        xmaX = np.max(data[:, 0])
        xmiX = np.min(data[:, 0])
        # for i in range(0, int(len(data[:, 0])/self.Nk1)):
            # ax.plot(data[(i)*self.Nk1:(i+1)*self.Nk1, 0], data[(i)*self.Nk1:(i+1)*self.Nk1, 1], color=(r, g, b), lw=1)
        # for i in range(0, len(data[:, 0])):
        #     ax.scatter(data[i, 0], data[i, 1], s=5, color=(rgbv[0, i], rgbv[1, i], rgbv[2, i]))
        cm = plt.cm.seismic
        # cm = plt.cm.jet
        ax.scatter(data[:, 0], data[:, 1], s=5, c=rgbv, cmap=cm)
        cb = plt.colorbar(mpl.cm.ScalarMappable(cmap=cm), ax=ax)
        ax.set_ylim(eng[0], eng[1])
        ax.set_xlim(xmiX, xmaX)
        ax.set_xticks([0, xmaX])
        ax.set_xticklabels(["$\Gamma$", "Z"])
        labels = ax.get_yticklabels() + ax.get_xticklabels()
        [label.set_fontname('Times New Roman') for label in labels]
        [label.set_fontsize(16) for label in labels]
        ax.axhline(y=0, xmin=0, xmax=100, linestyle='--', linewidth=2, color='cyan')
        plt.savefig(osp.join(self.file_path, "rib.png"), dpi=300, bbox_inches='tight')

def manipulate():
    figs = (10, 8)
    eng = [-1, 3]
    rib = ribbon()
    rib.plot_ribbon(figs, eng)

if __name__== "__main__":
    manipulate()

