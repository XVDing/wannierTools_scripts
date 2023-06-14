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
import func.funcs as F
import func.vb_bz as vb
curpath = os.path.abspath(os.path.dirname(__file__))
sys.path.append(curpath)

##################### process wannier tools ###########################
class wt():
    def __init__(self, filename, color_charity):
        self.soft_path = curpath
        self.file_path = os.getcwd()
        self.surface_card = F.read_surface_card(self.file_path)
        self.poscar = F.get_poscar(self.file_path, "POSCAR-rotated")
        self.font_size = 16
        self.filename = filename
        self.rotMat, self.ex, self.ey, self.ez = self.rot_3D()
        self.rot_pos, self.rep_surface, self.pos2D, self.rep_2D, self.n = self.matrix_data()
        self.realNodes = np.dot(F.get_nodes(self.file_path, filename), self.rotMat)
        self.arc_origionalP, self.arc_vector = F.read_arc_range(self.file_path, "WT.out")
        self.color_charity = color_charity
        self.pointsize = 10

    def norm(self, r):
        leng = [ii**2 for ii in r]
        rnew = np.sqrt(sum(leng))
        return rnew

    def visuzl_BZ(self):
        reciprocal_lattice_vectors=vb.read_poscar(self.file_path, 'POSCAR')   
        lattice_vectors=[np.array(reciprocal_lattice_vectors[0,:]), np.array(reciprocal_lattice_vectors[1,:]), np.array(reciprocal_lattice_vectors[2,:])]
        points, ridges, facets = vb.get_Wigner_Seitz_BZ(lattice_vectors)
        plt = vb.visualize_BZ_matplotlib(points, ridges, facets, reciprocal_lattice_vectors, self.file_path, self.filename, self.color_charity)
        return plt

    def rot_3D(self):
        pos = np.array(self.poscar).reshape(3, 3)
        urot = np.zeros(shape=(3, 3))
        Urot = []
        ex = pos[0]/self.norm(pos[0])
        Urot.append(ex)
        ez = F.chaCheng(pos[0], pos[1]) 
        ey = F.chaCheng(ez, ex)
        Urot.append(ey/self.norm(ey))
        Urot.append(ez/self.norm(ez))
        urot = np.array(Urot)
        return urot.T, ex, ey, ez

    def matrix_data(self):
        rot_mat = self.rotMat
        rot_pos = np.dot(np.array(self.poscar).reshape(3, 3), rot_mat)
        rep_surface_tmp = np.array(F.get_rec_surface(rot_pos.flatten())).reshape(3, 3)
        rep_surface = rep_surface_tmp[0:2]
        pos2D = np.around(np.array([rot_pos[0, 0:2], rot_pos[1, 0:2]]), 6)
        rep_2DTmp = F.get_rep2D(pos2D)[0:6]
        rep_2D = np.array(rep_2DTmp).reshape(2, 3)[:,0:2]
        ################ 法向量 ##################
        pos = np.array(self.poscar).reshape(3, 3)
        n_vertical = F.chaCheng(pos[0], pos[1])
        n = np.array(n_vertical / F.leng(n_vertical))
        return np.around(rot_pos, 6), np.around(rep_surface, 6), pos2D, rep_2D, n

    def write_node2d_ssband(self, nodes1):
        rep_surface_ni = np.linalg.inv(self.rep_surface[:, 0:2])
        nodes = np.dot(nodes1, rep_surface_ni)
        with open(osp.join(self.file_path, "node2d_SSband.dat"), 'w') as nd1:
            nd1.write("Nodes for 2D\n")
            for item in nodes:
                for i in range(2):
                    nd1.write(str(np.round(item[i], 6)).ljust(10, "0") + "   ")
                nd1.write("\n")
        nd1.close()

    def newcardirect(self):
        cards = []
        ######################## 旋转 ###############################
        nodes = F.get_nodes(self.file_path, self.filename)
        realnodes = np.dot(nodes, self.rotMat)
        ####################### 二维坐标系 ##########################
        n1 = np.array(F.chaCheng(self.rot_pos[0], self.rot_pos[1]))
        a1 = self.rot_pos[0]
        az = n1 
        a2 = np.array(F.chaCheng(az, a1))
        ####################### 二维坐标系 ##########################
        nodes_projected = []
        for jj in range(len(realnodes)):
            tmp = realnodes[jj] - F.dianCheng(realnodes[jj], n1)*n1
            nodes_projected.append(tmp)
        nodes_projected = np.array(nodes_projected)
        # nodes = np.transpose(np.dot(self.rotMat, np.transpose(nodes_projected)))
        # nodes = np.dot(nodes_projected, self.rotMat)
        
        ##################################################################
        for kk in range(len(realnodes)):
            ###################### 旋转后的坐标系下的分量 ###################################
            # x = F.dianCheng(self.rot_pos[0], nodes_projected[kk, :])/F.leng(self.rot_pos[0])
            # y = F.dianCheng(self.rot_pos[1], nodes_projected[kk, :])/F.leng(self.rot_pos[1])
            # z = F.dianCheng(self.rot_pos[2], nodes_projected[kk, :])/F.leng(self.rot_pos[2])
            ###################### 新坐标系下的分量 #########################################
            x = F.dianCheng(a1, nodes_projected[kk, :])/F.leng(a1)
            y = F.dianCheng(a2, nodes_projected[kk, :])/F.leng(a2)
            ################################################################################
            cards.append(np.array([x, y]))
        return np.array(cards)

    def write2D(self, filename):
        nodes2D = self.newcardirect()
        print("#################### 转换后的坐标 ########################")
        print(nodes2D)
        print("#########################################################")
        with open(osp.join(self.file_path, "nodes2D_arc.dat"), 'w') as arcnode:
            arcnode.write("Nodes for 2D arc"+"\n")
            for item1 in nodes2D:
                for i in range(2):
                    arcnode.write(str(np.round(item1[i], 6)).ljust(10, "0") + "   ")
                arcnode.write("\n")
        arcnode.close()
        ########################## nodes for ssband ###########################
        self.write_node2d_ssband(nodes2D)

    def stru_info(self):
        print("----------------------------------------------------------")
        print("")
        print("KPLANE_SLAB read from WT.out")
        print(self.arc_origionalP)
        print(self.arc_vector)
        print("")
        print("--- rotated matrix -----")
        print(self.rotMat)
        print(" ")
        print("--- rotated poscar (Z component of R1 and R1  are zero) -----")
        print(" ")
        print(np.around(self.rot_pos, 6))
        print(" ")
        print("--- rotated rep surface -----")
        print(" ")
        print(np.around(self.rep_surface, 6))
        print(" ")
        print("--- computed 2D rep  -----")
        print(" ")
        print(self.rep_2D)
        print(" ")
        print("--- read 2D rep from WT.out -----")
        print(" ")
        print(F.read_wtout(self.file_path))
        print(" ")
        print("---------------------------------------------------------")
        print("法向量: ")
        print(self.n)
        print("---------------------------------------------------------")
        print("直角坐标系下的 realnodes 坐标： ")
        print(self.realNodes)
        print("---------------------------------------------------------")

    def plot_arc(self, cm, filename="arc.dat_r", figname="arc_r.png"):
        x, y, z = F.data_pro(self.file_path, filename)
        plt.figure(figsize=(6, 6))
        ax = plt.subplot(111)
        # cm = plt.cm.seismic
        p = ax.scatter(x, y, s=5, c=z, cmap=cm, marker="o", vmin=np.min(z), vmax=np.max(z), alpha=0.8)
        cb = plt.colorbar(mpl.cm.ScalarMappable(cmap=cm), ax=ax)
        origional_points, arc_range = F.read_arc_range(self.file_path, "wt.in")
        #################################### nodes #############################
        ax.set_ylabel(r"$k_y$", fontsize=self.font_size, fontname="Times New Roman")
        ax.set_xlabel(r"$k_x$", fontsize=self.font_size, fontname="Times New Roman")
        nodes2d = np.loadtxt(osp.join(self.file_path, "nodes2D_arc.dat"), skiprows=1, dtype=np.float64)
        ##################### charity #########################################
        if osp.exists(osp.join(self.file_path, "wanniercenter3D_Weyl.dat")):
            charity = F.charity_data(self.file_path, "wanniercenter3D_Weyl.dat")
            for ii in range(len(charity)):
                if charity[ii] >= 1:
                    ax.text(nodes2d[ii, 0], nodes2d[ii, 1], str(ii), color=self.color_charity[0], fontsize=15, fontname="Times New Roman")
                    ax.scatter(nodes2d[ii, 0], nodes2d[ii, 1], s=self.pointsize, color=self.color_charity[0], marker="o")
                elif charity[ii] <= -1:
                    ax.text(nodes2d[ii, 0], nodes2d[ii, 1], str(ii), color=self.color_charity[1], fontsize=15, fontname="Times New Roman")
                    ax.scatter(nodes2d[ii, 0], nodes2d[ii, 1], s=self.pointsize, color=self.color_charity[1], marker="o")
        else:           
            for jj in range(len(nodes2d)):
                ax.scatter(nodes2d[jj, 0], nodes2d[jj, 1], s=self.pointsize, marker="o", color=self.color_charity[0])
                ax.text(nodes2d[jj, 0], nodes2d[jj, 1], str(jj), color=self.color_charity[0], fontsize=15, fontname="Times New Roman")
        ########################################################################
        # fig.colorbar(p, ax=ax)
        #ltick = 5
        #tick_x = np.array([item for item in np.arange(0, np.max(x), np.ceil(np.max(x)) / ltick)])
        #tick_x = np.append(tick_x, -tick_x).flatten()
        #tick_y = np.array([item for item in np.arange(0, np.max(y), np.ceil(np.max(y)) / ltick)])
        #tick_y = np.append(tick_y, -tick_y).flatten()
        #ax.set_xticks(tick_x)
        #ax.set_yticks(tick_y)
        labels = ax.get_yticklabels() + ax.get_xticklabels() + cb.ax.yaxis.get_ticklabels()
        [label.set_fontname('Times New Roman') for label in labels]
        [label.set_fontsize(16) for label in labels]
        ax.set_title(self.color_charity[0] + ": +1 or charity > 0")
        #vec = self.arc_origionalP[0]*self.rep_2D[0] + self.arc_origionalP[1]*self.rep_2D[1] + np.array([,])
        ax.set_xlim(np.min(x), np.max(x))
        ax.set_ylim(np.min(y), np.max(y))
        plt.savefig(figname, dpi=300, bbox_inches='tight')

def manipulate():
    # cm = plt.cm.jet
    # cm = plt.cm.gist_earth
    # cm = plt.cm.cubehelix
    # cm = plt.cm.seismic
    from matplotlib import pyplot as plt
    cm = plt.cm.viridis
    color_charity = ["red", "blue"]
    if os.path.exists("FBZ.dat"):
        filename_nodes="FBZ.dat"
    else:
        filename_nodes="Nodes.dat"
    filename_arc_r = "arc.dat_r"
    filename_arc_l = "arc.dat_l"
    wtt = wt(filename_nodes, color_charity)
    wtt.stru_info()
    wtt.write2D(filename_nodes) 
    wtt.plot_arc(cm, filename_arc_r, "arc_nodes_r.png")
    wtt.plot_arc(cm, filename_arc_l, "arc_nodes_l.png")
    plt = wtt.visuzl_BZ()
    #plt.show()

if __name__== "__main__":
    manipulate()

