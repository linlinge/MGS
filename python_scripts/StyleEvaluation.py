import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import open3d as o3d
import os
from sklearn.neighbors import KernelDensity
from numba import jit
from KDEpy import FFTKDE
from scipy.signal import find_peaks

# def OutlierGrade(input_path):
#     str_cmd="../build/MGS -i "+ input_path +" -o 1.csv -M OGEM"
#     os.system(str_cmd)
#     file=pd.read_csv("1.csv")
#     # og=np.asarray(file)
#     # print("Outlier Grade: %.4f" % og[0,0])
#     og=np.asarray(file)[:,0]
#     # print(np.mean(og))
#     x, y = FFTKDE(kernel="gaussian", bw=1).fit(og).evaluate()
#     # print(np.max(y))
#     print(len(y[y<np.max(y)*0.5])/len(x))
#     print(x[y==np.max(y)])
#     # peaks, _ = find_peaks(y, height=0.01)
#     # plt.plot(np.array([x[peaks]]), y[peaks]+0.3*np.ones([1,peaks.shape[0]]), "rv", markersize=8)    
#     plt.plot(x,y)
#     plt.savefig("1.png")
#     plt.legend(input_path.split("/")[-1])   
#     # plt.show()

def OutlierGrade(input_path,fig_name):
    str_cmd="../build/MGS -i "+ input_path +" -o 1.csv -M OutlierGrade"
    os.system(str_cmd)
    file=pd.read_csv("1.csv")
    og=np.asarray(file)[:,0]
    dty=np.asarray(file)[:,1]
    # dty=dty/np.sum(dty)
    # rst_og=og*dty
    # print("og min: %f"% np.min(og))
    # print("dty min: %f max: %f"% (np.min(dty),np.max(dty)))
    x, y = FFTKDE(kernel="gaussian", bw=1).fit(og).evaluate()
    # print(np.max(y))
    # print("og=%s" % str(len(y[y<np.max(y)*0.5])/len(x)))
    x_new=x[y>np.max(y)*0.5]
    og=(np.max(x_new)-np.min(x_new))/50
    print("og: %.2f" % (og))
    # print("og:  %.2f" % (x[y==np.max(y)]/np.max(x)))
    plt.plot(x,y)
    plt.savefig(fig_name)
    return og

def Homogeneity(input_path):
    str_cmd="../build/MGS -i "+ input_path +" -o 2.csv -M Homogeneity"
    os.system(str_cmd)
    file=pd.read_csv("2.csv")
    uty=np.asarray(file)[:,0]
    x, y = FFTKDE(kernel="gaussian", bw=0.01).fit(uty).evaluate()        
    plt.rcParams.update({'font.size': 22})
    plt.plot(x,y)
    peaks, _ = find_peaks(y, height=0.01)
    plt.plot(np.array([x[peaks]]), y[peaks]+0.3*np.ones([1,peaks.shape[0]]), "rv", markersize=8)  
    fname=input_path.split("/")[-1].split(".")[0]    
    plt.savefig("fig/Homogeneity_"+fname+".png")
    print("Hm: %d" % (len(peaks)))
    return len(peaks)

def Sigularity(input_path,fig_name="1.png"):
    str_cmd="../build/MGS -i "+ input_path +" -o 3.csv -M Singularity"
    os.system(str_cmd)
    df=pd.read_csv("3.csv")
    df=np.asarray(df)
    print("Sg: %.2f" % (df[0,0]))
    return df[0,0]
    # sg=df.to_numpy()[:,0]
    # x, y = FFTKDE(kernel="gaussian", bw=0.01).fit(sg).evaluate()
    # plt.plot(x,y)
    # plt.savefig(fig_name)
    # print(x[y==np.max(y)])
    # if sg[0,0] > 0.13:
    #     print("sg:  %.2f    --> 4th-style" % sg[0,0])
    # else:
    #     print("sg:  %.2f    --> 3rd-style" % sg[0,0])




if __name__ == '__main__':
    # input_path="/home/llg/Experiment/tanks_and_temples/Truck/Truck_COLMAP_clipped.ply"
    # rst_sg=Sigularity(input_path)
    # if rst_sg>4:
    #     print("4th-style")
    # else:
    #     og=OutlierGrade(input_path,"og.png")
    #     if og<0.1:
    #         print("1st-style")
    #     else:
    #         hm=Homogeneity(input_path)
    #         if hm<=1:
    #             print("2nd-style")
    #         else:
    #             print("3rd-style")

    # Style 1
    # Sigularity("/home/llg/dataset/3Dlib/Ignatius.ply")
    # Sigularity("/home/llg/dataset/3Dlib/Ignatius_COLMAP.ply")
    # OutlierGrade("/home/llg/dataset/3Dlib/Ignatius.ply")
    # Homogeneity("/home/llg/dataset/3Dlib/Ignatius.ply")

    # Style 2
    # Sigularity("/home/llg/dataset/3Dlib/Ignatius_COLMAP.ply")
    # OutlierGrade("/home/llg/dataset/3Dlib/Ignatius_COLMAP.ply")
    # Homogeneity("/home/llg/dataset/3Dlib/Ignatius_COLMAP.ply")

    # Style 3
    # Sigularity("/home/llg/dataset/3Dlib/torch_points.ply")
    # OutlierGrade("/home/llg/dataset/3Dlib/torch_points.ply")
    # Homogeneity("/home/llg/dataset/3Dlib/torch_points.ply")

    # Style 4
    # Sigularity("/home/llg/dataset/3Dlib/dog_colmap.ply")
    # OutlierGrade("/home/llg/dataset/3Dlib/dog_colmap.ply")
    # Homogeneity("/home/llg/dataset/3Dlib/dog_colmap.ply")

    # OutlierGrade("/home/llg/Human/H_FRM_0145_clipped.ply","1_1.png")
    # OutlierGrade("/home/llg/Human/fused_20200114_01_FRM_0073_1920.ply","1_2.png")
    # OutlierGrade("/home/llg/dataset/3Dlib/Ignatius_COLMAP.ply","1_2.png")
    # OutlierGrade("/home/llg/dataset/3Dlib/Ignatius_COLMAP_clipped.ply","1_2.png")
    # OutlierGrade("/home/llg/dataset/3Dlib/cat_colmap.ply","1_2.png")
    # OutlierGrade("/home/llg/dataset/3Dlib/Ignatius.ply","1_2.png")
    OutlierGrade("/home/llg/dataset/3Dlib/Ignatius_COLMAP.ply","1_2.png")
    # Homogeneity("/home/llg/dataset/3Dlib/Ignatius.ply")
    # Homogeneity("/home/llg/dataset/3Dlib/Ignatius_COLMAP.ply")
    # Homogeneity("/home/llg/dataset/3Dlib/Barn_COLMAP.ply")
    # Homogeneity("/home/llg/dataset/3Dlib/torch_points.ply")
    # Sigularity("/home/llg/dataset/3Dlib/Meetingroom_COLMAP.ply")


    # Sigularity("/home/llg/dataset/3Dlib/cat_colmap.ply","3_1.png")
    # Sigularity("/home/llg/dataset/3Dlib/dog_colmap.ply","3_2.png")
    # Sigularity("/home/llg/dataset/3Dlib/Ignatius.ply","3_3.png")
    # Sigularity("/home/llg/dataset/3Dlib/Ignatius_COLMAP.ply","3_4.png")