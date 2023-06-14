
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
character = ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P", "Q", "R", "S", "T", 
            "U", "V", "W", "X", "Y", "Z", "a", "b", "c", "d", "e", "h", "i", "j", "k", "l", "m", "n", "o", "p", "q", "r", "s", "t", "u", "v", "w", "x", "y", "z"]
def gen_kp_ssband(filepath, filename="node2d_SSband.dat"):
    file = osp.join(filepath, filename)
    file1 = osp.join(filepath, "kpath.dat")
    data = np.loadtxt(file, skiprows=1, dtype=float)
    len_data = len(data[:, 0])
    if len_data <= 52:
        with open(file1, "w") as kpt:
            kpt.write(str(len_data-1) + "\n")
            for i in range(len(data[:, 0])-1):
                kpt.write(str(character[i]) + " " + str(np.around(data[i, 0], 8)).ljust(10, " ") + "  " + 
                    str(np.around(data[i, 1], 8)).ljust(10, " ") + "   " + character[i+1] + " " + str(np.around(data[i+1, 0], 8)).ljust(10, " ") 
                    + "  " + str(np.around(data[i+1, 1], 8)).ljust(10, " ") + "\n")
        kpt.close()

def get_ticks(filepath, filename="surfdos_l.gnu"):
    data = open(osp.join(filepath, filename), "r")
    data_lines = data.readlines()
    for i in range(len(data_lines)):
        if "yrange" in data_lines[i]:
            if ")" not in data_lines[i]:
                data_ll = data_lines[i+1] + data_lines[i+2]
                ticks = str(str(data_ll.split("(")[1]).split(")")[0]).split(",")
                ticks_label = [eval(ticks[i].split()[0]) for i in range(len(ticks))]
                ticks_data = [np.around(float(ticks[i].split()[1]), 8) for i in range(len(ticks))]
                print(ticks_label)
                print(ticks_data)
                break
            else:
                ticks = str(str(data_lines[i+1].split("(")[1]).split(")")[0]).split(",")
                ticks_label = [eval(ticks[i].split()[0]) for i in range(len(ticks))]
                ticks_data = [np.around(float(ticks[i].split()[1]), 8) for i in range(len(ticks))]
                print(ticks_label)
                print(ticks_data)
                break
    return np.array(ticks_label), ticks_data
############################## data processing ##############################################
def data_pro(filepath="./", filename=""):
    # data = np.loadtxt(filename, skiprows=0, dtype=np.float64)
    # X = data[:, 0]
    # Y = data[:, 1]
    # z = data[:, 2]
    nedos = 0
    data = np.loadtxt(filename, skiprows=0, dtype=np.float64)
    for i in range(0, len(data[:, 0])):
        if data[i, 0] != 0:
            nedos = i
            break
    x = data[0:-1:nedos, 0]
    y = data[:nedos, 1]
    z = np.transpose((data[:, 2].reshape(len(x), len(y))))
    X, Y = np.meshgrid(x, y)
    return X, Y, z

##################################### plot ##########################################
def plot_arc(filepath="./", filename="dos.dat_r", figname="ss_band_r.png", energy_range=[], font_size=16):
    x, y, z = data_pro(filepath, filename)
    plt.figure(figsize=(8, 6))
    ax = plt.subplot(111)
    cm = plt.cm.jet
    # cm = plt.cm.gist_earth
    # cm = plt.cm.cubehelix
    # cm = plt.cm.seismic
    # cm = plt.cm.seismic
    # cm = plt.cm.viridis
    # cm = plt.cm.hot
    print(np.min(z))
    print(np.max(z))
    p = ax.pcolormesh(x, y, z, cmap=cm, vmin=0, vmax=10, shading='flat', antialiased=True)
    # p = ax.pcolormesh(x, y, z, cmap=cm, vmin=0, vmax=4, shading='flat', antialiased=True)
    # p = ax.scatter(x, y, s=5, c=z, cmap=cm, marker="o", vmin=np.min(z), vmax=np.max(z), alpha=0.8)
    cb = plt.colorbar(p)
    ########################################################################
    ax.set_ylim(energy_range)
    ticks_label, ticks_data = get_ticks(filepath, "surfdos_l.gnu")
    ticks_data[-1] = np.max(x)
    for i in range(0,len(ticks_label)):
        ax.axvline(ticks_data[i],linewidth=0.8,linestyle="--", color='gray')

    for i in range(len(ticks_label)):
        if ticks_label[i] == "G":
            ticks_label[i] = u"Γ"
    ax.set_xticklabels(ticks_label, fontsize=font_size, fontdict=font)
    # ax.set_xticks([])
    # ax.set_xticklabels([])
    # ax.axhline(y=eng, xmin=0, xmax=100, linestyle='--', linewidth=1, color='white')
    ax.set_xlim(np.min(x), np.max(x))
    # ax.set_ylim(energy_range[0], energy_range[1])
    # ax.set_title("Weyl Point 1 and 2")
    labels = ax.get_yticklabels() + ax.get_xticklabels() + cb.ax.yaxis.get_ticklabels()
    # labels = ax.get_yticklabels() + cb.ax.yaxis.get_ticklabels()
    [label.set_fontproperties(font_properties) for label in labels]
    [label.set_fontsize(24) for label in labels]
    plt.savefig(figname, dpi=300, bbox_inches='tight')
    plt.close()

if __name__== "__main__":
    file_path = os.getcwd()
    energy_range = [56, 61]
    print("---------------------------------------------------------")
    print("begin plot ssband !")
    plot_arc(filepath=file_path, filename="dos.dat_r", figname="ss_band_r.png", energy_range=energy_range, font_size=16)
    plot_arc(filepath=file_path, filename="dos.dat_l", figname="ss_band_l.png", energy_range=energy_range, font_size=16)


