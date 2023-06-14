
import matplotlib.pyplot as plt
import numpy as np
import os.path as osp
import matplotlib as mpl
import os, sys
curpath = os.path.abspath(os.path.dirname(__file__))
sys.path.append(curpath)

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

def gen_kp_ssband_biaoji(filepath, filename="node2d_SSband.dat", biaoji=[]):
    file = osp.join(filepath, filename)
    file1 = osp.join(filepath, "kpath.dat")
    tmp = np.loadtxt(file, skiprows=1, dtype=float)
    data_tmp = []
    for ii in range(len(biaoji)):
        data_tmp.append(tmp[biaoji[ii], :])
    data = np.array(data_tmp)
    len_data = len(data)
    if len_data <= 52:
        with open(file1, "w") as kpt:
            kpt.write(str(len_data-1) + "\n")
            for i in range(len(data[:, 0])-1):
                kpt.write(str(character[i]) + " " + str(np.around(data[i, 0], 8)).ljust(10, " ") + "  " + 
                    str(np.around(data[i, 1], 8)).ljust(10, " ") + "   " + character[i+1] + " " + str(np.around(data[i+1, 0], 8)).ljust(10, " ") 
                    + "  " + str(np.around(data[i+1, 1], 8)).ljust(10, " ") + "\n")
        kpt.close()
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
    plt.figure(figsize=(6, 6))
    ax = plt.subplot(111)
    # cm = plt.cm.jet
    # cm = plt.cm.gist_earth
    # cm = plt.cm.cubehelix
    # cm = plt.cm.seismic
    # cm = plt.cm.seismic
    cm = plt.cm.viridis
    p = ax.pcolormesh(x, y, z, cmap=cm, vmin=-3, vmax=6, shading='flat', antialiased=True)
    # p = ax.scatter(x, y, s=5, c=z, cmap=cm, marker="o", vmin=np.min(z), vmax=np.max(z), alpha=0.8)
    cb = plt.colorbar(p)
    ########################################################################
    ax.set_xlim(np.min(x), np.max(x))
    ax.set_ylim(np.min(y), np.max(y))
    # ax.set_ylim(energy_range[0], energy_range[1])
    labels = ax.get_yticklabels() + ax.get_xticklabels() + cb.ax.yaxis.get_ticklabels()
    [label.set_fontname('Times New Roman') for label in labels]
    [label.set_fontsize(16) for label in labels]
    plt.savefig(figname, dpi=800, bbox_inches='tight')

if __name__== "__main__":
    file_path = os.getcwd()
    fname = "node2d_SSband.dat"
    energy_range = [-0.5, 0.5]
    print("---------------------------------------")
    print("Select a numbers: ")
    print("----   1. gkp     2. plot ssband:  ----")
    print("---------------------------------------")
    print("---------------------------------------------------------")
    num = int(input())
    print("")
    print("--------------------------------------")
    if num == 1:
        print("---------------------------------------")
        print("Select a numbers: ")
        print("---      1. defaule                ---")
        print("---      2. give biaoji index      ---")
        print("---------------------------------------")
        seclect = int(input())
        if seclect == 2:
            print("begin generate kpt !")
            biaoji_tmp = input("give biaoji index, separate by blank: ").split()
            biaoji = [int(kk) for kk in biaoji_tmp]
            gen_kp_ssband_biaoji(file_path, fname, biaoji)
        else:
            print("begin generate kpt !")
            gen_kp_ssband(file_path, fname)
    else:
        print("begin plot ssband !")
        cmdd = "cp " + osp.join(curpath, "ssband.sh") + " ."
        os.system(cmdd)
        plot_arc(filepath=file_path, filename="dos.dat_r", figname="ss_band_r.png", energy_range=energy_range, font_size=16)
        plot_arc(filepath=file_path, filename="dos.dat_l", figname="ss_band_l.png", energy_range=energy_range, font_size=16)


