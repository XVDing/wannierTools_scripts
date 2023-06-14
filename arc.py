
import matplotlib.pyplot as plt
import numpy as np
import os.path as osp
import os, sys
curpath = os.path.abspath(os.path.dirname(__file__))
sys.path.append(curpath)
import matplotlib as mpl
from matplotlib import font_manager

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
def data_pro(filepath="./", filename=""):
    data = np.loadtxt(filename, skiprows=0, dtype=np.float64)
    X = data[:, 0]
    Y = data[:, 1]
    z = data[:, 2]
    # nedos = 0
    # data = np.loadtxt(filename, skiprows=0, dtype=np.float64)
    # for i in range(0, len(data[:, 0])):
    #     if data[i, 0] != data[0, 0]:
    #         nedos = i
    #         break
    # x = data[0:-1:nedos, 0]
    # y = data[:nedos, 1]
    # z = np.transpose((data[:, 2].reshape(len(x), len(y))))
    # X, Y = np.meshgrid(x, y)
    return X, Y, z

##################################### plot ##########################################
def plot_arc(filepath="./", filename="dos.dat_r", figname="ss_band_r.png", font_size=16):
    x, y, z = data_pro(filepath, filename)
    plt.figure(figsize=(8, 6))
    ax = plt.subplot(111)
    # cm = plt.cm.jet
    # cm = plt.cm.gist_earth
    # cm = plt.cm.cubehelix
    # cm = plt.cm.seismic
    # cm = plt.cm.seismic
    # cm = plt.cm.viridis
    cm = plt.cm.twilight_shifted
    # cm = plt.cm.hot
    # p = ax.pcolormesh(x, y, z, cmap=cm, vmin=np.min(z)+3, vmax=np.max(z)+4, shading='flat', antialiased=True)
    # p = ax.pcolormesh(x, y, z, cmap=cm, vmin=0, vmax=4, shading='flat', antialiased=True)
    p = ax.scatter(x, y, s=5, c=z, cmap=cm, marker="o", vmin=np.min(z), vmax=np.max(z), alpha=0.8)
    cb = plt.colorbar(p)
    ########################################################################
    ax.set_ylim(np.min(y), np.max(y))
    ax.set_xlim(np.min(x), np.max(x))
    ax.set_xticks([])
    ax.set_yticks([])
    plt.savefig(figname, dpi=600, bbox_inches='tight')
    plt.close()

if __name__== "__main__":
    file_path = os.getcwd()
    print("---------------------------------------------------------")
    print("begin plot ssband !")
    plot_arc(filepath=file_path, filename="arc.dat_r", figname="arc-r.png", font_size=16)
    plot_arc(filepath=file_path, filename="arc.dat_l", figname="arc-l.png", font_size=16)


