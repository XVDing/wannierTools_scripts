import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import os.path as osp
import matplotlib as mpl
import os, sys
import func.funcs as F
from matplotlib import font_manager

curpath = os.path.abspath(os.path.dirname(__file__))
sys.path.append(curpath)
# ------------------ font setup ----------------------#
fontpath='/work/wangr/dxy/scripts/font/arial.ttf'
font_properties = font_manager.FontProperties(fname=fontpath)
styles = ['normal', 'italic', 'oblique']
weights = ['light', 'normal', 'medium', 'semibold', 'bold', 'heavy', 'black']
font = {'fontproperties': font_properties,  # 'family': 'Liberation Sans'
		'style': styles[0],
		'color': 'black',
		'weight': weights[1],
		'size': 18, 
        }
font_small = {'fontproperties': font_properties,  # 'family': 'Liberation Sans'
		'style': styles[0],
		'color': 'black',
		'weight': weights[1],
		'size': 12, 
        }
plt.rc('font',family=font, size=14)
matplotlib.use('agg')

def plot_bulkBand(filepath, filename, eng):
    Nk = F.get_NK(filepath, "WT.out")
    file = osp.join(filepath, filename)
    data = np.loadtxt(file, skiprows=0, dtype=float)
    ax = plt.subplot(111)
    x = data[:Nk[0], 0]
    # print(len(data[:, 1])/Nk[0])
    energy = data[:, 1].reshape(int(len(data[:, 1])/Nk[0]), Nk[0])
    for i in range(len(energy[:, 0])):
        ax.plot(x, energy[i, :], color="blue")
    ticks_label, ticks_data = F.get_ticks(filepath, "bulkek.gnu")
    ticks_data[-1] = np.max(x)
    for i in range(0,len(ticks_label)):
        ax.axvline(ticks_data[i],linewidth=0.8,linestyle="--", color='gray')
    ax.set_xticks(ticks_data)
    ax.set_xticklabels(ticks_label, fontsize=16)
    ax.set_xlim(np.min(x), np.max(x))
    ax.set_ylim(eng)
    ax.set_ylabel("Energy (eV)")
    plt.savefig(osp.join(filepath, "bulkband.png"), dpi=300, bbox_inches='tight')
    plt.close()

def get_index(energy, yuzhi):
    index = []
    for i in range(len(energy[:, 0])):
        if np.abs(np.max(energy[:, i]) -0) < yuzhi and np.abs(np.min(energy[:, i+1])-0) < 0.1:
            index.append([i, i+1])
    return index

def bulk_band_parity(filepath, filename, figss, linesize, linesize_kuang, eng, font_size, yuzhi=0.1):
    Nk = F.get_NK(filepath, "WT.out")
    file = osp.join(filepath, filename)
    data = np.loadtxt(file, skiprows=0, dtype=float)
    ########################### parity_rb ##########################
    plt.rcParams['xtick.direction'] = 'in'#将x周的刻度线方向设置向内
    plt.rcParams['ytick.direction'] = 'in'#将y轴的刻度方向设置向内
    plt.figure(figsize=(figss))
    ax = plt.subplot(111)
    ax.spines['bottom'].set_linewidth(str(linesize_kuang))
    ax.spines['top'].set_linewidth(str(linesize_kuang))
    ax.spines['left'].set_linewidth(str(linesize_kuang))
    ax.spines['right'].set_linewidth(str(linesize_kuang))
    x = data[:Nk[0], 0]
    # print(len(data[:, 1])/Nk[0])
    energy = data[:, 1].reshape(int(len(data[:, 1])/Nk[0]), Nk[0])
    index_vbm_cbm = get_index(energy, yuzhi)
    for i in range(len(energy[:, 0])):
        ax.plot(x, energy[i, :], linewidth=linesize, color="blue")
    ticks_label, ticks_data = F.get_ticks(filepath, "bulkek.gnu")
    ticks_data[-1] = np.max(x)
    for i in range(0,len(ticks_label)):
        ax.axvline(ticks_data[i],linewidth=0.8,linestyle="--", color='gray')
    # ax.set_xticks(ticks_data)
    # ax.set_xticklabels(ticks_label, fontsize=16)
    ax.set_xticks([])
    ytick = [eng[0] + i*0.2 for i in range(3)]
    ax.set_yticks(ytick)
    ax.set_xlim(np.min(x)+np.max(x)*2/5, np.max(x)*2.8/5)
    # ax.set_xlim(np.min(x), np.max(x))
    ax.set_yticks([])
    ax.set_ylim(eng)
    plt.ylim(eng)
    # ax.set_ylabel("Energy (eV)")
    label_x = ax.get_xticklabels()
    [label.set_fontproperties(font_properties) for label in label_x]
    [label.set_fontsize(font_size) for label in label_x]
    label_y = ax.get_yticklabels() 
    [label.set_fontproperties(font_properties) for label in label_y]
    [label.set_fontsize(font_size-2) for label in label_y]
    plt.savefig(osp.join(filepath, "band.png"), dpi=300, bbox_inches='tight')
    plt.close()

if __name__== "__main__":
    file_path = os.getcwd()
    figs=(3, 3.5)
    energy_range = [0, 64]
    ls = 1.0
    ls_kuang = 0.5
    font_s = 14
    padinch = 0.01
    yuzhi_weyl = 0.1
    plot_bulkBand(file_path, "bulkek.dat", energy_range)
    bulk_band_parity(file_path, "bulkek.dat", figs, ls, ls_kuang, energy_range, font_s, yuzhi_weyl)
 