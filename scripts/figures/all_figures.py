from barplot import main as barplot_main
from boxplot import main as boxplot_main
from broken_regions import main as broken_regions_main
from heatmap import main as heatmap_main
from sub_families import main as sub_families_main
from venn_diagram import main as venn_diagram_main


path = "/home/jaimy/output/human/bcr_v2/"

barplot_main(path)
boxplot_main(path)
broken_regions_main(path)
sub_families_main(path)
venn_diagram_main(path, "bcr")
heatmap_main(path)