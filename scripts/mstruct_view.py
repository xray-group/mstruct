#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import math
import sys


def mstruct_view(filename, Ihkl_file_list=[], xl=[]):
    d = np.loadtxt(filename + '.dat')
    x = d[:,0]
    y = d[:,1]
    dy = d[:,1]
    
    if len(xl) == 0:
        xl.append(x.min())
        xl.append(x.max())
    
    ind = np.where((x >= xl[0]) & (x <= xl[1]))[0]
    
    x = x[ind]
    y = y[ind]
    dy = dy[ind]
    
    #label_str = 'Intensity (counts)'
    label_str = 'Intensity (arbitrary units)'
    
    fig = plt.figure()
    fig.set_size_inches(18/2.54, 14/2.54, forward=True)
    ax = fig.add_subplot(111)
    ax.set_xlabel('$2\Theta (^\circ)$')
    ax.set_ylabel(label_str)
    ax.set_xlim([xl[0], xl[1]])
    ax.set_title(filename)
    
    ax.plot(x, np.sqrt(d[ind,1]), color='C3', label='data')
    ax.plot(x, np.sqrt(d[ind,2]), color='C0', label='fit')
    ax.plot(x, np.sign(d[ind,3]) * np.sqrt(np.abs(d[ind,3])), color='C2', label='data-fit')
    
    ax.legend(loc="upper right", frameon=False)
    
    if len(Ihkl_file_list) > 0:
        #colors = ['darkred', 'darkgreen', 'blue', 'black']
        colors = ['C0', 'C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'C7']
        
        for i, iname in enumerate(Ihkl_file_list):
            hkl = []
            with open('Ihkl_' + iname + '.dat', 'r') as f:
                for line in f:
                    if line in ['\n', '\r\n']:
                        break
                    if line[0] == '#':
                        continue
                    vals = [float(s) for s in line.split()] 
                    th2 = vals[3]
                    Fhkl = vals[4]
                    if (th2 >= xl[0]) and (th2 <= xl[1]) and (Fhkl > 1e-7):
                        txt = '- {:d} {:d} {:d}'.format(int(vals[0]), int(vals[1]), int(vals[2]))
                        #print txt
                        yhkl = math.sqrt(1.2 * np.interp(th2, d[:,0], d[:,1]));
                        ax.text(th2, yhkl, txt, fontsize=9, color=colors[i], horizontalalignment='center', verticalalignment='bottom', rotation='vertical')
    
    
    plt.tight_layout(pad=0.0)
    plt.show()


if __name__ == "__main__":
    filename = sys.argv[1]
    Ihkl_file_list = []
    xl = []
    if len(sys.argv) > 2:
        Ihkl_file_list = [s for s in sys.argv[2].replace('[','').replace(']','').split(',')]
        if len(sys.argv) > 3:
            xl = [int(s) for s in sys.argv[3].replace('[','').replace(']','').split(',')]
    
    mstruct_view(filename, Ihkl_file_list, xl)
