import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import open3d as o3d
import os
from sklearn.neighbors import KernelDensity
from numba import jit
from KDEpy import FFTKDE
from scipy.signal import find_peaks

def OutlierGrade(input_path):
    str_cmd="../build/MGS -i "+ input_path +" -o 1.csv -M OGEM"
    os.system(str_cmd)
    file=pd.read_csv("1.csv")
    og=np.asarray(file)
    print("Outlier Grade:   "+str(og[0,0]))

def Uniformity(input_path):
    str_cmd="../build/MGS -i "+ input_path +" -o 2.csv -M UEM"
    os.system(str_cmd)
    file=pd.read_csv("2.csv")
    uty=np.asarray(file)[:,0]
    x, y = FFTKDE(kernel="gaussian", bw=0.01).fit(uty).evaluate()
    
    plt.rcParams.update({'font.size': 22})
    plt.plot(x,y)
    peaks, _ = find_peaks(y, height=0.01)
    plt.plot(np.array([x[peaks]]), y[peaks]+0.3*np.ones([1,peaks.shape[0]]), "rv", markersize=8)    
    plt.show()

def Sigularity(input_path):
    str_cmd="../build/MGS -i "+ input_path +" -o 3.csv -M SEM"
    os.system(str_cmd)
    file=pd.read_csv("3.csv")
    sg=np.asarray(file)
    if sg[0,0] > 0.1:
        print("Singularity: %.4f    --> Recommand Style 4" % sg[0,0])
    else:
        print("Singularity: %.4f    --> Recommand Style 3" % sg[0,0])


# OutlierGrade("/home/llg/dataset/3Dlib/cat_colmap.ply")
# Uniformity("/home/llg/dataset/3Dlib/cat_colmap.ply")
Sigularity("/home/llg/dataset/3Dlib/Meetingroom_COLMAP.ply")
# Sigularity("/home/llg/dataset/3Dlib/cat_colmap.ply")

