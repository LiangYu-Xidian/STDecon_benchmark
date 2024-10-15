import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pylab import mpl
mpl.rcParams['font.sans-serif'] = ['STZhongsong']    # default
mpl.rcParams['axes.unicode_minus'] = False
plt.figure(figsize=(4, 4))

#100_jsd
x =["CARD","Cell2location","CellDART","DestVI","DSTG","novoSpaRc","RCTD","SD2","Seurat","SpatialDecon","SpatialDWLS","SPOTlight","SpSeudoMap","STdeconvolve","stereoscope","Tangram"]
y =[0.22609,0.13689,0.30077,0.08687,0.27066,0.231,0.1557,0.27664,0.38973,0.23057,0.54305,0.17663,0.42232,0.16546,0.1452,0.21364]

colors = ['#FA2F25', '#EB9885', '#F7A403', '#6CC0D4', '#94FF08', '#487FDE', '#00FFE8', '#9AA028', '#021DFF', '#D500FF', '#8B0000', '#FF00CA', '#FF8C00', '#8B008B', '#00CED1', '#FFD700']
plt.barh(x,width=y,align="center",color=colors)
# plt.title("JSD↓",loc="center",fontsize=12)
plt.xlabel('JSD',fontsize=12)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.savefig('jsd_100.svg',bbox_inches='tight')

#100_rmse
x =["CARD","Cell2location","CellDART","DestVI","DSTG","novoSpaRc","RCTD","SD2","Seurat","SpatialDecon","SpatialDWLS","SPOTlight","SpSeudoMap","STdeconvolve","stereoscope","Tangram"]
y =[0.188018969,0.142529568,0.244112824,0.130482858,0.217187879,0.200589515,0.183944525,0.233247865,0.281246246,0.222259648,0.354754906,0.169089953,0.243042658,0.163399785,0.174482772,0.190193936]

colors = ['#FA2F25', '#EB9885', '#F7A403', '#6CC0D4', '#94FF08', '#487FDE', '#00FFE8', '#9AA028', '#021DFF', '#D500FF', '#8B0000', '#FF00CA', '#FF8C00', '#8B008B', '#00CED1', '#FFD700']
plt.barh(x,width=y,align="center",color=colors)
#plt.title("RMSE↓",loc="center",fontsize=12)
plt.xlabel('RMSE',fontsize=12)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.savefig('rmse_100.svg',bbox_inches='tight')

#100_pcc
x =["CARD","Cell2location","CellDART","DestVI","DSTG","novoSpaRc","RCTD","SD2","Seurat","SpatialDecon","SpatialDWLS","SPOTlight","SpSeudoMap","STdeconvolve","stereoscope","Tangram"]
y =[0.402589534,0.735581403,0.290421787,0.760797836,0.162844111,0.434121494,0.517507334,0.195503863,0.054657519,0.381842249,0.131199658,0.54498295,0.187651152,0.699064058,0.603596692,0.377687249]
colors = ['#FA2F25', '#EB9885', '#F7A403', '#6CC0D4', '#94FF08', '#487FDE', '#00FFE8', '#9AA028', '#021DFF', '#D500FF', '#8B0000', '#FF00CA', '#FF8C00', '#8B008B', '#00CED1', '#FFD700']
plt.barh(x,width=y,align="center",color=colors)
#plt.title("PCC↑",loc="center",fontsize=12)
plt.xlabel('PCC',fontsize=12)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.savefig('pcc_100.svg',bbox_inches='tight')

#100_spearman
x =["CARD","Cell2location","CellDART","DestVI","DSTG","novoSpaRc","RCTD","SD2","Seurat","SpatialDecon","SpatialDWLS","SPOTlight","SpSeudoMap","STdeconvolve","stereoscope","Tangram"]
y =[0.286385898,0.697655053,0.305818113,0.786999315,0.114321702,0.32387467,0.603099077,0.256087049,0.064055379,0.44013788,0.113522745,0.472935935,0.009058852,0.5680785,0.627521344,0.407266889]
colors = ['#FA2F25', '#EB9885', '#F7A403', '#6CC0D4', '#94FF08', '#487FDE', '#00FFE8', '#9AA028', '#021DFF', '#D500FF', '#8B0000', '#FF00CA', '#FF8C00', '#8B008B', '#00CED1', '#FFD700']
plt.barh(x,width=y,align="center",color=colors)
#plt.title("Spearman↑",loc="center",fontsize=12)
plt.xlabel('Spearman',fontsize=12)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.savefig('spearman_100.svg',bbox_inches='tight')