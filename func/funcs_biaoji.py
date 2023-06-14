# -*- coding: utf-8 -*-
"""
Created on 19:38 07-11 2022

@author: XY Ding
mail to: dxy_vasp@163.com
python3: arc.py
"""
import math
import numpy as np
from matplotlib import pyplot as plt
import matplotlib as mpl
import os, sys 
import os.path as osp

# get reciprocal unit cell
def chaCheng(a1, a2):
    c = []
    c1 = float(a1[1] * a2[2] - a1[2] * a2[1])
    c.append(c1)
    c2 = float(a1[2] * a2[0] - a1[0] * a2[2])
    c.append(c2)
    c3 = float(a1[0] * a2[1] - a1[1] * a2[0])
    c.append(c3)
    return c

def dianCheng(b1, b2):
    d1 = float(b1[0] * b2[0])
    d2 = float(b1[1] * b2[1])
    d3 = float(b1[2] * b2[2])
    d = d1 + d2 + d3
    return d
def get_poscar(filepath="", filename='POSCAR'):
    file = osp.join(filepath, filename)
    poscar = open(file, 'r')
    poscar_lines = poscar.readlines()
    pos = []
    for flag_lines in range(2, 5):
        for i in range(0, 3):
            pos.append(float(poscar_lines[flag_lines].split()[i]))
    return pos

def get_rec3D(filepath, filename):
    pos = get_poscar(filepath, filename)
    rep = []
    a1 = [pos[0], pos[1], pos[2]]
    a2 = [pos[3], pos[4], pos[5]]
    a3 = [pos[6], pos[7], pos[8]]
    volume = dianCheng(a1, chaCheng(a2, a3))     # volume is a1 . (a2 X a3)
    scalar = 2 * math.pi / volume
    for i in range(0, 3):
        b1 = scalar*chaCheng(a2, a3)[i]
        rep.append(b1)
    for j in range(0, 3):
        b2 = scalar*chaCheng(a3, a1)[j]
        rep.append(b2)
    for k in range(0, 3):
        b3 = scalar * chaCheng(a1, a2)[k]
        rep.append(b3)
    return rep

def get_rec(pos):
    rep = []
    a1 = [pos[0], pos[1], pos[2]]
    a2 = [pos[3], pos[4], pos[5]]
    a3 = [pos[6], pos[7], pos[8]]
    volume = dianCheng(a1, chaCheng(a2, a3))     # volume is a1 . (a2 X a3)
    scalar = 2 * math.pi / volume
    for i in range(0, 3):
        b1 = scalar*chaCheng(a2, a3)[i]
        rep.append(b1)
    for j in range(0, 3):
        b2 = scalar*chaCheng(a3, a1)[j]
        rep.append(b2)
    for k in range(0, 3):
        b3 = scalar * chaCheng(a1, a2)[k]
        rep.append(b3)
    return rep

def get_rec_surface(pos):
    rep = []
    a1 = [pos[0], pos[1], pos[2]]
    a2 = [pos[3], pos[4], pos[5]]
    a3 = [0, 0, 1]
    volume = dianCheng(a1, chaCheng(a2, a3))     # volume is a1 . (a2 X a3)
    scalar = 2 * math.pi / volume
    for i in range(0, 3):
        b1 = scalar*chaCheng(a2, a3)[i]
        rep.append(b1)
    for j in range(0, 3):
        b2 = scalar*chaCheng(a3, a1)[j]
        rep.append(b2)
    for k in range(0, 3):
        b3 = scalar * chaCheng(a1, a2)[k]
        rep.append(b3)
    return rep

def dianCheng2D(a1, a2):
    d1 = float(a1[0] * a2[0])
    d2 = float(a1[1] * a2[1])
    d = d1 + d2
    return d
def get_pos(filepath, filename, surface_card):
    pos = np.array(get_poscar(filepath, filename)).reshape(3, 3)
    pos_surface_card = np.dot(surface_card, pos)
    return pos_surface_card
def leng(a):
    temp = 0
    for i in range(3):
        temp = temp + a[i]**2
    return np.sqrt(float(temp))
def leng_2d(a):
    temp = 0
    for i in range(2):
        temp = temp + a[i]**2
def get_rep_surf(filepath, filename, surface_card):
    pos_surf_card = get_pos(filepath, filename, surface_card)
    a1 = pos_surf_card[0]
    a2 = pos_surf_card[1]
    a3_temp = np.array(chaCheng(a1, a2))
    a3 = a3_temp/leng(a3_temp)
    scal = 2 * math.pi / dianCheng(a1, chaCheng(a2, a3))
    b1 = scal * np.array(chaCheng(a2, a3))
    b2 = scal * np.array(chaCheng(a3, a1))
    b = np.array([b1, b2])
    return b

def get_rep2D(pos):
    rep = []
    a1 = [pos[0][0], pos[0][1], 0]
    a2 = [pos[1][0], pos[1][1], 0]
    a3 = [0, 0, 1]
    volume = dianCheng(a1, chaCheng(a2, a3))
    scalar = 2 * math.pi / volume
    for i in range(0, 3):
        b1 = scalar*chaCheng(a2, a3)[i]
        rep.append(b1)
    for j in range(0, 3):
        b2 = scalar*chaCheng(a3, a1)[j]
        rep.append(b2)
    for k in range(0, 3):
        b3 = scalar * chaCheng(a1, a2)[k]
        rep.append(b3)
    return rep

def leng(a):
    temp = 0
    for i in range(3):
        temp = temp + a[i]**2
    return np.sqrt(float(temp))

def leng_2d(a):
    temp = 0
    for i in range(2):
        temp = temp + a[i]**2
    return np.sqrt(float(temp))

def get_theta(a, b):
    def jia_jiao(a, b):
        dotab = np.dot(a, b)
        la = leng(a)
        lb = leng(b)
        theta = np.arccos(float(dotab / (la * lb)))
        return theta
    def jia_jiao_2d(a, b):
        dotab = np.dot(a, b)
        la = leng_2d(a)
        lb = leng_2d(b)
        theta = np.arccos(float(dotab / (la * lb)))
        return theta
    theta = 0
    if len(a) == 3:
        theta = jia_jiao(a, b)
    else:
        theta = jia_jiao_2d(a,b)
    return theta

def read_wtout(filepath):
    file = osp.join(filepath, "WT.out")
    wt = open(file, 'r').readlines()
    k = []
    for i in range(len(wt)):
        if "Ka2" in wt[i]:
            k.append(wt[i+1].split())
            k.append(wt[i+2].split())
    return np.array(k, dtype=np.float64)

def rot_matrix(theta, a):
    cos = np.cos(theta)
    sin = np.sin(theta)
    rot = [[cos + (1-cos)*a[0]**2, (1-cos)*a[0]*a[1]-sin*a[2],
            (1-cos)*a[0]*a[2] + sin * a[1]],
           [(1-cos)*a[0]*a[1] + sin*a[2], cos + (1-cos)*a[1]**2,
            (1-cos)*a[1]*a[2] - sin * a[0]],
           [(1-cos)*a[0]*a[2]-sin*a[1], (1-cos)*a[1]*a[2]+sin*a[0],
            cos + (1-cos)*a[2]**2]]
    return rot

def read_surface_card(filepath):
    file = osp.join(filepath, "WT.out")
    wt = open(file, 'r').readlines()
    surface = []
    for i in range(len(wt)):
        if "1st vector on surface" in wt[i]:
            surface.append(wt[i].split(":")[1].split())
            surface.append(wt[i+1].split(":")[1].split())
    return np.array(surface, dtype=np.float64)

def data_pro(filepath, filename="arc.dat_r"):
    file = osp.join(filepath, filename)
    data = np.loadtxt(file, skiprows=0, dtype=np.float64)
    X = data[:, 0]
    Y = data[:, 1]
    z = data[:, 2]
    return X, Y, z
def get_nodes(filepath, filename):
    file = osp.join(filepath, filename)
    node = np.loadtxt(file, skiprows=2, dtype=float)
    nodes = node[:, 0:3]
    return nodes

def charity_data(filepath, filename="wanniercenter3D_Weyl.dat"):
    file = osp.join(filepath, filename)
    data = open(file, 'r').readlines()
    charity = data[0].split()[2:]
    if osp.exists(osp.join(filepath, "FBZ.dat")):
        index = []
        data_nodes = np.loadtxt(osp.join(filepath, "Nodes.dat"), skiprows=2, dtype=float)
        data_FBZ = np.loadtxt(osp.join(filepath, "FBZ.dat"), skiprows=2, dtype=float)
        for j in range(len(data_nodes[:, 0])):
            for i in range(len(data_FBZ[:, 0])):
                if data_nodes[j, 5] == data_FBZ[i, 5] and data_nodes[j, 6] == data_FBZ[i, 6] and data_nodes[j, 7] == data_FBZ[i, 7]:
                    index.append(j)
                    break
        charity_FBZ = []
        for ii in range(len(index)):
            charity_FBZ.append(charity[ii])
        return np.array(charity_FBZ, dtype=float)
    else:
        return np.array(charity, dtype=float)
def get_nodes_biaoji(filepath, filename, biaoji):
    file = osp.join(filepath, filename)
    node = np.loadtxt(file, skiprows=2, dtype=float)
    nodes = []
    for item in biaoji:
        nodes.append(node[item, 0:3])
    return np.array(nodes)

def charity_biaoji(filepath="./", filename="wanniercenter3D_Weyl.dat", biaoji=[]):
    charity = charity_data(filepath, filename="wanniercenter3D_Weyl.dat")
    charity_biaoji = []
    for ii in biaoji:
        charity_biaoji.append(charity[ii])
    return np.array(charity_biaoji)

def read_arc_range(filepath, filename="WT.out"):
    file = osp.join(filepath, filename)
    data = open(file, 'r').readlines()
    origional_points = []
    first_vector = []
    second_vector = []
    for ii in range(len(data)):
        if 'K2D_start' in data[ii]:
            origional_points = [float(item) for item in str(data[ii].split(":")[1]).split()]
            first_vector = [float(item) for item in str(data[ii+1].split(":")[1]).split()]
            second_vector = [float(item) for item in str(data[ii+2].split(":")[1]).split()]
            break
    return np.array(origional_points), np.array([first_vector, second_vector])

character = ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P", "Q", "R", "S", "T", "U", "V", "W", "X", "Y", "Z", "a", "b", "c", "d", "e", "f"]
def gen_kp_ssband(filepath):
    file = osp.join(filepath, "node2d_SSband.dat")
    file1 = osp.join(filepath, "kpath.dat")
    data = np.loadtxt(file, skiprows=1, dtype=float)
    len_data = len(data[:, 0])
    if len_data <= 26:
        kpt = open(file1, "w")
        kpt.write(str(len_data-1) + "\n")
        for i in range(len(data[:, 0])-1):
            kpt.write(str(character[i]) + " " + str(np.around(data[i, 0], 4)).ljust(6, " ") + "  " + 
            str(np.around(data[i, 1], 4)).ljust(6, " ") + "   " + character[i+1] + " " + str(np.around(data[i+1, 0], 4)).ljust(6, " ") 
                    + "  " + str(np.around(data[i+1, 1], 4)).ljust(6, " ") + "\n")
        kpt.close()
############################## data processing ##############################################
def data_pro(filepath="./", filename=""):
    file = osp.join(filepath, filename)
    data = np.loadtxt(file, skiprows=0, dtype=np.float64)
    X = data[:, 0]
    Y = data[:, 1]
    z = data[:, 2]
    return X, Y, z

##################################### plot ##########################################
def plot_ssband(filepath="./", filename="arc.dat_r", figname="ss_band_r.png", font_size=16):
    x, y, z = data_pro(filepath, filename)
    plt.figure(figsize=(6, 6))
    ax = plt.subplot(111)
    cm = plt.cm.seismic
    p = ax.scatter(x, y, s=5, c=z, cmap=cm, marker="o", vmin=np.min(z), vmax=np.max(z), alpha=0.8)
    cb = plt.colorbar(mpl.cm.ScalarMappable(cmap=cm), ax=ax)
    ########################################################################
    ax.set_xlim(np.min(x), np.max(x))
    ax.set_ylim(np.min(y), np.max(y))
    labels = ax.get_yticklabels() + ax.get_xticklabels() + cb.ax.yaxis.get_ticklabels()
    [label.set_fontname('Times New Roman') for label in labels]
    [label.set_fontsize(16) for label in labels]
    plt.savefig(figname, dpi=800)

if __name__== "__main__":
    # surf_card = [[-1, 0, 1], [0, 1, 0]]
    #print(rot)
    if os.path.exists("FBZ.dat"):
        file_name="FBZ.dat"
    else:
        file_name="Nodes.dat"
    print("test")

