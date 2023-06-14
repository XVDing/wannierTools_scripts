import numpy as np, numpy.linalg as npl
from scipy.spatial import Voronoi
from itertools import product
from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d import proj3d
import matplotlib as mpl
import os.path as osp
import math
mpl.rcParams['font.size'] = 16.
##################################################### basic function ###################################################
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

def get_reciprocal(filepath="", filename='POSCAR'):
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
    return np.array(rep).reshape(3, 3)

def read_high_symmetry_points(filepath, filename = "KPOINTS"):
    file = osp.join(filepath, filename)
    high_kpt_point = []
    kpt = []
    kpoints = open(file, 'r').read().strip('\n').splitlines()
    ''' 
    filter() 函数用于过滤序列，过滤掉不符合条件的元素，返回由符合条件元素组成的新列表。
    该接收两个参数，第一个为函数，第二个为序列，序列的每个元素作为参数传递给函数进行判断，
    然后返回 True 或 False，最后将返回 True 的元素放到新列表中。
    '''
    kpoints = list(filter(str.strip, kpoints))
    tmp = ''
    for flag_hk in range(4, len(kpoints)):
        kpt_hp_data = kpoints[flag_hk].split('!')[1].strip()
        if kpoints[flag_hk].split('!')[0] == tmp:
            continue
        else:
            kpt_hpd_data = list(filter(str.strip, kpoints[flag_hk].split('!')[0].split(' ')))
            if kpt_hp_data == "\Gamma":
                kpt_hp_data = u"Γ"
            high_kpt_point.append(kpt_hp_data)
            kpt.append(kpt_hpd_data)
        tmp = kpoints[flag_hk].split('!')[0]
    kpt = np.array(kpt, dtype=float)
    for i in range(0, len(high_kpt_point)-1):
        if high_kpt_point[i] is high_kpt_point[i + 1]:
            high_kpt_point[i] = ' '
    high_kpt_point = list(filter(str.strip, high_kpt_point))
    print(kpt)
    print(high_kpt_point)
    return kpt, high_kpt_point
########################################################################################################################
def read_poscar(filepath, filename="POSCAR", species=None):
    file = osp.join(filepath, filename)
    poscar = open(file,'r')
    title = poscar.readline().strip()
    scale = float(poscar.readline().strip())
    s = float(scale)
    lattice_vectors = [[ float(v) for v in poscar.readline().split() ],
            [ float(v) for v in poscar.readline().split() ],
            [ float(v) for v in poscar.readline().split() ]]
    lattice_vectors = np.array(lattice_vectors)
    reciprocal_lattice_vectors= np.linalg.inv(lattice_vectors).T
    reciprocal_lattice_vectors=reciprocal_lattice_vectors*np.pi*2
    return reciprocal_lattice_vectors

def is_greek_alphabets(klabels):
    Greek_alphabets=['Alpha','Beta','Gamma','Delta','Epsilon','Zeta','Eta','Theta', 'Iota','Kappa','Lambda','Mu','Nu','Xi','Omicron','Pi','Rho','Sigma','Tau','Upsilon','Phi','Chi','Psi','Pega']
    group_labels=[]
    for i in range(len(klabels)): 
        klabel=klabels[i] 
        for j in range(len(Greek_alphabets)):
            if (klabel.find(Greek_alphabets[j].upper())>=0):
                latex_exp=r''+'$\\'+str(Greek_alphabets[j])+'$'
                klabel=klabel.replace(str(Greek_alphabets[j].upper()),str(latex_exp))
        if (klabel.find('_')>0):
           n=klabel.find('_')
           klabel=klabel[:n]+'$'+klabel[n:n+2]+'$'+klabel[n+2:]
        group_labels.append(klabel)
    klabels=group_labels
    return klabels

def read_kpath(filepath, filename="KPOINTS"):
    file = osp.join(filepath, filename)
    kpath=np.loadtxt(file, dtype=np.string_,skiprows=4)
    #print(kpath)
    kpath_labels = kpath[:,3].tolist()
    kpath_labels = [i.decode('utf-8','ignore') for i in kpath_labels]
    for i in range(len(kpath_labels)):
           if kpath_labels[i]=="Gamma":
                   kpath_labels[i]=u"Γ"
    kpaths=np.zeros((len(kpath_labels),3),dtype=float)
    for i in range(len(kpath_labels)):
        kpaths[i,:]=\
        [float(x) for x in kpath[i][0:3]]
    return kpath_labels, kpaths

def get_Wigner_Seitz_BZ(lattice_vectors):
# Inspired by http://www.thp.uni-koeln.de/trebst/Lectures/SolidState-2016/wigner_seitz_3d.py 
# Inspired by https://github.com/QijingZheng/VASP_FermiSurface/blob/master/fs.py 
    latt = []
    prefactors = [0., -1., 1.]
    for p in prefactors:
        for u in lattice_vectors:
            latt.append(p * u)
    lattice = []
    for vs in product(latt, latt, latt):
        a = vs[0] + vs[1] + vs[2]
        if not any((a == x).all() for x in lattice):
            lattice.append(a)
    voronoi = Voronoi(lattice)
    bz_facets = []
    bz_ridges = []
    bz_vertices = []
    for pid, rid in zip(voronoi.ridge_points, voronoi.ridge_vertices):
        if(pid[0] == 0 or pid[1] == 0):
            bz_ridges.append(voronoi.vertices[np.r_[rid, [rid[0]]]])
            bz_facets.append(voronoi.vertices[rid])
            bz_vertices += rid
    bz_vertices = list(set(bz_vertices))
    return voronoi.vertices[bz_vertices], bz_ridges, bz_facets

class Arrow3D(FancyArrowPatch):
    def __init__(self, xs, ys, zs, *args, **kwargs):
        FancyArrowPatch.__init__(self, (0, 0), (0, 0), *args, **kwargs)
        self._verts3d = xs, ys, zs
    def draw(self, renderer):
        xs3d, ys3d, zs3d = self._verts3d
        xs, ys, zs = proj3d.proj_transform(xs3d, ys3d, zs3d, self.axes.M)
        self.set_positions((xs[0], ys[0]), (xs[1], ys[1]))
        FancyArrowPatch.draw(self, renderer)

def data_process_nodes(filepath, filename="Nodes.dat"):
    file = osp.join(filepath, filename)
    node = np.loadtxt(file, skiprows=2, dtype=float)
    kx = node[:, 0]
    ky = node[:, 1]
    kz = node[:, 2]
    return kx, ky, kz

def charity_data(filepath, filename="wanniercenter3D_Weyl.dat"):
    file = osp.join(filepath, filename)
    data = open(file, 'r').readlines()
    charity = data[0].split()[2:]
    if osp.exists("FBZ.dat"):
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

def charity_biaoji(filepath="./", filename="wanniercenter3D_Weyl.dat", biaoji=[]):
    charity = charity_data(filepath, filename)
    charity_biaoji = []
    for itt in biaoji:
        charity_biaoji.append(charity[itt])
    return np.array(charity_biaoji)

def visualize_BZ_matplotlib(points,ridges,facets,reciprocal_lattice_vectors, filepath, filename, biaoji, color_charity):
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection
    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(111, projection='3d', proj_type='ortho')
    basis_vector_clrs = ['r', 'g', 'b']
    basis_vector_labs = ['b$_1$', 'b$_2$', '$b_3$']
    for ii in range(3):
        arrow = Arrow3D([0,reciprocal_lattice_vectors[ii, 0]], [0,reciprocal_lattice_vectors[ii, 1]], [0,reciprocal_lattice_vectors[ii, 2]],
                color=basis_vector_clrs[ii], mutation_scale=10,lw=1,arrowstyle="->")
        ax.add_artist(arrow)
        ax.text(reciprocal_lattice_vectors[ii, 0], reciprocal_lattice_vectors[ii, 1],reciprocal_lattice_vectors[ii, 2],
                basis_vector_labs[ii], fontsize=16, fontname="Times New Roman")
        for ir in ridges:
            ax.plot(ir[:, 0], ir[:, 1], ir[:, 2], color='k', lw=1.5,alpha=0.5)
    k1, k2, k3 = data_process_nodes(filepath, filename)
    #ax.scatter(k1, k2, k3, c='r', marker='o',s=10,alpha=0.8)
    if osp.exists(osp.join(filepath, "wanniercenter3D_Weyl.dat")):
        charity = charity_biaoji(filepath, "wanniercenter3D_Weyl.dat", biaoji)
        for jj in range(len(charity)):
            if charity[jj] <= -1:
                ax.text(k1[biaoji[jj]], k2[biaoji[jj]], k3[biaoji[jj]], str(biaoji[jj]), color=color_charity[1], fontsize=15, fontname="Times New Roman")
                ax.scatter(k1[biaoji[jj]], k2[biaoji[jj]], k3[biaoji[jj]], marker='o', color=color_charity[1],s=15,alpha=0.8)
            elif charity[jj] >= 1:
                ax.text(k1[biaoji[jj]], k2[biaoji[jj]], k3[biaoji[jj]], str(biaoji[jj]), color=color_charity[0],  fontsize=15, fontname="Times New Roman")
                ax.scatter(k1[biaoji[jj]], k2[biaoji[jj]], k3[biaoji[jj]], marker='o', color=color_charity[0],s=15,alpha=0.8)
    else:
        for ii in range(len(k1)):
            ax.scatter(k1[biaoji[jj]], k2[biaoji[jj]], k3[biaoji[jj]], c=color_charity[0], marker='o',s=15,alpha=0.8)
            ax.text(k1[biaoji[jj]], k2[biaoji[jj]], k3[biaoji[jj]], str(ii), color=color_charity[0],  fontsize=15, fontname="Times New Roman")
    ax.set_axis_off()
    ax.view_init(elev=16, azim=-70)
    plt.savefig('bz_biaoji.png', dpi=300, bbox_inches='tight')
    ax.view_init(elev=0, azim=45)
    plt.savefig('bz_side_biaoji.png', dpi=300, bbox_inches='tight')
    ax.view_init(elev=90, azim=-45)
    plt.savefig('bz_top_biaoji.png', dpi=300, bbox_inches='tight')
    return plt
    
if __name__ == "__main__":   
   reciprocal_lattice_vectors=read_poscar('POSCAR')   
   lattice_vectors=[np.array(reciprocal_lattice_vectors[0,:]),np.array(reciprocal_lattice_vectors[1,:]),np.array(reciprocal_lattice_vectors[2,:])]
   #print(lattice_vectors)
   points, ridges, facets = get_Wigner_Seitz_BZ(lattice_vectors)
   visualize_BZ_matplotlib(points, ridges, facets, reciprocal_lattice_vectors)
