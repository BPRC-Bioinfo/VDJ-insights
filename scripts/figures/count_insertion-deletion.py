import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd


path = "/home/jaimy/output/human/bcr_v2/"
data_novel = pd.read_excel(f"{path}/annotation/annotation_report_novel.xlsx")
data_novel = data_novel[data_novel['Segment'].isin(['V', 'D', 'J'])]
grouped_data = data_novel.groupby(['Region', 'Segment', 'Short name', 'BTOP']).size().reset_index(name='Count').sort_values(by='Count', ascending=False)

print(grouped_data)
count_samples = data_novel['Sample'].nunique()

plt.figure(figsize=(12, 6))
plot = sns.catplot(
    data=grouped_data,
    x="Segment",
    y="Count",
    col="Region",
    col_wrap=2,
    height=5,
    aspect=1.5,
    palette="Set2",
    #kind="violin",
    #sharex=False,
)


plot.set_titles("{col_name}", size=12)
plot.fig.suptitle(f"Regional patterns of mutation frequency across samples (N={count_samples})", fontsize=16)
plot.set_axis_labels("Segment", "Count", fontsize=14)
plt.tight_layout()
plt.savefig(f"{path}/insertion_deletion_count_novel.png", dpi=600)