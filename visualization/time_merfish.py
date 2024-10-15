import pandas as pd
import matplotlib.pyplot as plt

# read data
df = pd.read_csv(r"D:\MERFISH\time.csv")
colors = ['#FA2F25', '#EB9885', '#F7A403', '#6CC0D4', '#94FF08', '#487FDE', '#00FFE8', '#9AA028', '#021DFF', '#D500FF', '#8B0000', '#FF00CA', '#FF8C00', '#8B008B', '#00CED1', '#FFD700']
# setting
methods = df['Unnamed: 0'].tolist()  # x label
time_points = df.columns[1:].tolist()
performance_data = df.iloc[:, 1:].values  # performance


fig, axs = plt.subplots(1, 1, figsize=(5, 4.5))

# for each method
for i in range(len(methods)):
    axs.plot(time_points, performance_data[i], marker='o',linewidth=1, markersize=5, color=colors[i])


axs.set_xlabel('bin size /μm', fontsize=12)
axs.set_ylabel('Computation Time /s',fontsize=12)
axs.tick_params(axis='both', which='major', labelsize=12)
axs.tick_params(axis='both', which='minor', labelsize=6)
#for spine in axs.spines.values():
#spine.set_linewidth(3)

# show
plt.tight_layout()

plt.savefig("time_MERFISH.svg",dpi=1000,bbox_inches='tight')
#plt.plot(df.columns[1:], marker='o', linewidth=2.5, markersize=10, color=colors)
#plt.show()
