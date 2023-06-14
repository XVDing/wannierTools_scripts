# -*- coding: utf-8 -*-
"""
Created on 19:38 07-11 2022

@author: XY Ding
mail to: dxy_vasp@163.com
python3: arc.py
"""
import numpy as np
from matplotlib import pyplot as plt
import matplotlib as mpl
import os, sys
import os.path as osp
from matplotlib import font_manager
# fontpath='/work/wangr/dxy/scripts/font/times.ttf'
fontpath='/work/wangr/dxy/scripts/font/arial.ttf'
# fontpath="D:\\jianguoyun\\Scripts\\font\\times.ttf"
# fontpath="D:\\jianguoyun\\Scripts\\font\\arial.ttf"
# ------------------ font setup ----------------------#
font_properties = font_manager.FontProperties(fname=fontpath)
styles = ['normal', 'italic', 'oblique']
weights = ['light', 'normal', 'medium', 'semibold', 'bold', 'heavy', 'black']
font = {'fontproperties': font_properties,  # 'family': 'Liberation Sans'
		'style': styles[0],
		'color': 'black',
		'weight': weights[1],
		'size': 16, 
        }
plt.rc('font',family=font, size=20)

#################################################################################
def get_arc_data(filepath, filename="WT.out"):
    file = osp.join(filepath, filename)
    data = open(file, 'r').readlines()
    kplane_slab = {}
    vector = []
    for i in range(len(data)):
        if 'K2D_start' in data[i]:
            kplane_slab['orp'] = np.array([float(item) for item in str(data[i].split(':')[1]).split()])
            vector.append(np.array([float(item1) for item1 in str(data[i+1].split(':')[1]).split()]))
            vector.append(np.array([float(item2) for item2 in str(data[i+2].split(':')[1]).split()]))
            kplane_slab['vector'] = vector
            return kplane_slab

def read_wtout_Ka2(filepath):
    file = osp.join(filepath, "WT.out")
    wt = open(file, 'r').readlines()
    k = []
    for i in range(len(wt)):
        if "Ka2" in wt[i]:
            k.append(wt[i+1].split())
            k.append(wt[i+2].split())
    return np.array(k, dtype=np.float64)

def trans_nodes_arc(filepath, filename='nodes2D_arc.dat', filename2="nodes2D_arc_transformed.dat"):
    kplane_card = get_arc_data(os.getcwd())
    ka2_kb2 = read_wtout_Ka2(os.getcwd())
    vector_direct = np.dot(kplane_card['vector'], ka2_kb2)
    origional_point = np.dot(kplane_card['orp'], vector_direct)
    x_max = origional_point[0] + vector_direct[0, 0]
    y_max = origional_point[1] + vector_direct[1, 1]
    max_point = [x_max, y_max]
    print('********************* kplane card dirt *************************')
    print(origional_point)
    print(vector_direct)
    print("xrange:  " )
    print(origional_point[0], x_max)
    print("yrange: ")
    print(origional_point[1], y_max)
    print('********************* kplane card dirt *************************')
    file = osp.join(filepath, filename)
    data_nodes_arc = np.loadtxt(file, skiprows=1, dtype=float)
    save_nodes = open(osp.join(filepath, filename2), 'w')
    save_nodes.write('Nodes for 2D arc')
    save_nodes.write('\n')
    for i in range(len(data_nodes_arc)):
    # for i in range(3):
        print('*********** The ' + str(i) + ' Point ! *************')
        back = False 
        tmp_nodes = data_nodes_arc[i, :]
        while not back:
            if tmp_nodes[0] > origional_point[0] and tmp_nodes[0] < x_max:
                if tmp_nodes[1] > origional_point[1] and tmp_nodes[1] < y_max:
                    for k in range(len(tmp_nodes)):
                        save_nodes.write(str(np.round(tmp_nodes[k], 12)).rjust(15) + "  ")
                    save_nodes.write('\n')  
                    back = True 
                    # print("good")
                else:
                    max_multiply = max([abs(int((tmp_nodes[1] - origional_point[1])/(y_max-origional_point[1]))), abs(int((tmp_nodes[1] - y_max)/(y_max-origional_point[1])))])
                    if max_multiply == abs(int((tmp_nodes[1] - origional_point[1])/(y_max-origional_point[1]))):
                        mul = int((tmp_nodes[1] - origional_point[1])/(y_max-origional_point[1]))
                    else:
                        mul = int((tmp_nodes[1] - y_max)/(y_max-origional_point[1]))
                    tmp_nodes = tmp_nodes - mul * vector_direct[1, :]
                    print("********* 1 *******")
                    print(tmp_nodes)
                    # break 
            elif tmp_nodes[0] < origional_point[0] or tmp_nodes[0] > x_max:
                if tmp_nodes[1] > origional_point[1] and tmp_nodes[1] < y_max:
                    max_multiply1 = max([abs(int((tmp_nodes[0] - origional_point[0])/(x_max-origional_point[0]))), abs(int((tmp_nodes[0] - x_max)/(x_max-origional_point[0])))])
                    if max_multiply1 == abs(int((tmp_nodes[0] - origional_point[0])/(x_max-origional_point[0]))):
                        mul1 = int((tmp_nodes[0] - origional_point[0])/(x_max-origional_point[0]))
                    else:
                        mul1 = int((tmp_nodes[0] - x_max)/(x_max-origional_point[0]))
                    tmp_nodes = tmp_nodes - mul1 * vector_direct[0, :] 
                    print("********* 2 *******")
                    print(tmp_nodes)
            else:
                max_multiplyx = max([abs(int((tmp_nodes[0] - origional_point[0])/(x_max-origional_point[0]))), abs(int((tmp_nodes[0] - x_max)/(x_max-origional_point[0])))])
                if max_multiplyx == abs(int((tmp_nodes[0] - origional_point[0])/(x_max-origional_point[0]))):
                    mulx = int((tmp_nodes[0] - origional_point[0])/(x_max-origional_point[0]))
                else:
                    mulx = int((tmp_nodes[0] - x_max)/(x_max-origional_point[0]))
                max_multiplyy = max([abs(int((tmp_nodes[1] - origional_point[1])/(y_max-origional_point[1]))), abs(int((tmp_nodes[1] - y_max)/(y_max-origional_point[1])))])
                if max_multiplyy == abs(int((tmp_nodes[1] - origional_point[1])/(y_max-origional_point[1]))):
                    muly = int((tmp_nodes[1] - origional_point[1])/(y_max-origional_point[1]))
                else:
                    muly = int((tmp_nodes[1] - y_max)/(y_max-origional_point[1]))
                tmp_nodes = tmp_nodes - mulx * vector_direct[0, :] - muly * vector_direct[1, :] 
                print("********* 3 *******") 
                print(tmp_nodes)
                # break 
    print('*********** End transform Points ! *************')
    save_nodes.close() 
    os.system('mv nodes2D_arc_transformed.dat nodes2D_arc.dat')     

def plot_arc(cm, file_path, filename="arc.dat_r", figname="arc_r.png", font_size=18, color_charity=['red', 'blue']):
    import func.funcs as F
    x, y, z = F.data_pro(file_path, filename)
    plt.rcParams['xtick.direction'] = 'in' #将x周的刻度线方向设置向内
    plt.rcParams['ytick.direction'] = 'in' #将y轴的刻度方向设置向内
    plt.figure(figsize=(4.5, 4.0))
    ax = plt.subplot(111)
    # cm = plt.cm.seismic
    p = ax.scatter(x, y, s=5, c=z, cmap=cm, marker="o", vmin=np.min(z)+2, vmax=np.max(z)+3, alpha=0.8)
    # cb = plt.colorbar(mpl.cm.ScalarMappable(cmap=cm), ax=ax)
    # origional_points, arc_range = F.read_arc_range(file_path, "wt.in")
    #################################### nodes #############################
    # ax.set_ylabel(r"$k_y$", fontsize=font_size, fontdict=font")
    # ax.set_xlabel(r"$k_x$", fontsize=font_size, fontname="Times New Roman")
    nodes2d = np.loadtxt(osp.join(file_path, "nodes2D_arc.dat"), skiprows=1, dtype=np.float64)
    ##################### charity #########################################
    ps = 8
    if osp.exists(osp.join(file_path, "wanniercenter3D_Weyl.dat")):
        charity = F.charity_data(file_path, "wanniercenter3D_Weyl.dat")
        for ii in range(len(charity)):
            if charity[ii] >= 1:
                # ax.text(nodes2d[ii, 0], nodes2d[ii, 1], str(ii), color=color_charity[0], fontsize=15, fontname="Times New Roman")
                ax.scatter(nodes2d[ii, 0], nodes2d[ii, 1], s=ps, color=color_charity[0], marker="o")
            elif charity[ii] <= -1:
                # ax.text(nodes2d[ii, 0], nodes2d[ii, 1], str(ii), color=color_charity[1], fontsize=15, fontname="Times New Roman")
                ax.scatter(nodes2d[ii, 0], nodes2d[ii, 1], s=ps, color=color_charity[1], marker="o")
    else:           
        for jj in range(len(nodes2d)):
            ax.scatter(nodes2d[jj, 0], nodes2d[jj, 1], s=ps, marker="o", color=color_charity[0])
            # ax.text(nodes2d[jj, 0], nodes2d[jj, 1], str(jj), color=color_charity[0], fontsize=15, fontname="Times New Roman")
    ########################################################################
    # fig.colorbar(p, ax=ax)
    #ltick = 5
    #tick_x = np.array([item for item in np.arange(0, np.max(x), np.ceil(np.max(x)) / ltick)])
    #tick_x = np.append(tick_x, -tick_x).flatten()
    #tick_y = np.array([item for item in np.arange(0, np.max(y), np.ceil(np.max(y)) / ltick)])
    #tick_y = np.append(tick_y, -tick_y).flatten()
    #ax.set_xticks(tick_x)
    #ax.set_yticks(tick_y)
    # labels = ax.get_yticklabels() + ax.get_xticklabels() + cb.ax.yaxis.get_ticklabels()
    labels = ax.get_yticklabels() + ax.get_xticklabels() 
    [label.set_fontproperties(font_properties) for label in labels]
    # [label.set_fontsize(font_size) for label in labels]
    # ax.set_title(color_charity[0] + ": +1 or parity > 0")
    #vec = self.arc_origionalP[0]*self.rep_2D[0] + self.arc_origionalP[1]*self.rep_2D[1] + np.array([,])
    ax.set_xlim(np.min(x), np.max(x))
    ax.set_ylim(np.min(y), np.max(y))
    plt.savefig(figname, dpi=600, bbox_inches='tight')
    plt.close()
def manipulate():
    # cm = plt.cm.jet
    # cm = plt.cm.gist_earth
    # cm = plt.cm.cubehelix
    # cm = plt.cm.seismic
    from matplotlib import pyplot as plt
    cm = plt.cm.viridis
    # cm = plt.cm.hot
    color_charity = ["red", "Fuchsia"]
    filename_arc_r = "arc.dat_r"
    filename_arc_l = "arc.dat_l"
    fonts = 16
    plot_arc(cm, os.getcwd(), filename_arc_r, "trans_arc_nodes_r.png", fonts, color_charity)
    plot_arc(cm, os.getcwd(), filename_arc_l, "trans_arc_nodes_l.png", fonts, color_charity)

if __name__== "__main__":
    filepath = os.getcwd()
    trans_nodes_arc(filepath)           ############# do you need transform ? 
    print("********** Begin plotting ! ***********")
    manipulate()
    print("*********** End plotting ! ***********")


