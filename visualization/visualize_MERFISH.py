import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf


fig, axs = plt.subplots(2, 1, figsize=(11, 4.7))

# read_data
datasets = [
    ("D:\MERFISH\MERFISH_jsd.csv", "Cell Type Prediction (JSD)"),
    ("D:\MERFISH\MERFISH_rmse.csv", "Cell Type Prediction (RMSE)"),
]
colors = ['#FA2F25', '#EB9885', '#F7A403', '#6CC0D4', '#94FF08', '#487FDE', '#00FFE8', '#9AA028', '#021DFF', '#D500FF', '#8B0000', '#FF00CA', '#FF8C00', '#8B008B', '#00CED1', '#FFD700']
merged_data = []
for file_path, title in datasets:
    df = pd.read_csv(file_path)
    merged_data.append((title, df))


for i, (title, df) in enumerate(merged_data):
    ax = axs[i]
    for index, row in df.iterrows():
        color = colors[index % len(colors)]
        ax.plot(df.columns[1:], row[1:], marker='o', label=row[0], linewidth=1, markersize=5, color=color)
    # ax.set_xlabel('bin size /μm', fontsize=16, fontweight="bold")
    if i == 0:
        ax.set_ylabel('JSD', fontsize=12)
    elif i == 1:
        ax.set_ylabel('RMSE', fontsize=12)
    # elif i == 2:
    #     ax.set_ylabel('Computation Time /s', fontsize=16, fontweight="bold")
    ax.tick_params(axis='both', which='major', labelsize=12)
    ax.tick_params(axis='both', which='minor', labelsize=12)


    #ax.text(-0.03, 1.02, chr(65 + i), transform=ax.transAxes,
    #        fontsize=16, fontweight='bold', va='top', ha='right')



for ax in axs:
    ax.margins(x=0.01)


plt.xlabel('bin size /μm', fontsize=12)

plt.savefig("MERFISH_robutness.svg", bbox_inches='tight', dpi=1000)