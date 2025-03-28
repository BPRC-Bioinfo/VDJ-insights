
class NetworkBuilder:
    def __init__(self, alignment, tooltip_data, left_flanking_gene, right_flanking_gene, accessions, color_theme):
        self.alignment = alignment
        self.tooltip_data = tooltip_data
        self.left_flanking_gene = left_flanking_gene
        self.right_flanking_gene = right_flanking_gene
        self.accessions = accessions
        self.nodes_dict = {}
        self.links = []
        self.node_ids = set()
        self.y_positions, self.x_step, self.x_offset, self.node_radius, self.segment_color = self.configure_visualization(color_theme)


    def configure_visualization(self, color_theme):
        y_positions = {"haplotype1": 250, "common": 300, "haplotype2": 350}
        x_step = 35
        x_offset = 100
        node_radius = 10
        color_themes = {
            "default": {"V": "#ffdfba", "D": "#baffc9", "J": "#bae1ff", "C": "#7a6fac", "N/A": "#ffffff"},
            "easy_eye": {"V": "#E69F00", "D": "#56B4E9", "J": "#009E73", "C": "#F5F5DC", "N/A": "#ffffff"},
            "tol": {"V": "#DDCC77", "D": "#332288", "J": "#88CCEE", "C": "#F5F5DC", "N/A": "#ffffff"},
            "IBM": {"V": "#785EF0", "D": "#DC267F", "J": "#FE6100", "C": "#F5F5DC", "N/A": "#ffffff"}
        }
        return y_positions, x_step, x_offset, node_radius, color_themes.get(color_theme, "default")


    def build_network(self):
        left_gene_x = self.x_offset - self.x_step
        self.add_node(self.left_flanking_gene, 0, "white", left_gene_x, self.y_positions["common"], sample=self.accessions[0])
        hap1_prev = hap2_prev = self.nodes_dict[0]["index"]

        index = 0
        for idx, (seg1, seg2) in enumerate(self.alignment):
            if seg1 is not None:
                seg1 = seg1.split(" (")[0]
            if seg2 is not None:
                seg2 = seg2.split(" (")[0]

            x = idx * self.x_step + self.x_offset
            if seg1 == seg2 and seg1 is not None:
                index = self.add_common_node(index, seg1, x, hap1_prev, hap2_prev)
                hap1_prev = hap2_prev = index
            else:
                if seg1 is not None:
                    index = self.add_haplotype_node(index, seg1, x, "haplotype1", hap1_prev)
                    hap1_prev = index
                if seg2 is not None:
                    index = self.add_haplotype_node(index, seg2, x, "haplotype2", hap2_prev)
                    hap2_prev = index

        x = len(self.alignment) * self.x_step + self.x_offset
        self.add_node(self.right_flanking_gene, index + 1, "white", x, self.y_positions["common"], sample=self.accessions[0])
        self.add_link(hap1_prev, index + 1)
        self.add_link(hap2_prev, index + 1)

        return {
            "title": "version 8",
            "nodes": list(self.nodes_dict.values()),
            "links": self.links,
            "x_step": self.x_step,
            "x_offset": self.x_offset,
            "y_positions": self.y_positions,
            "node_radius": self.node_radius,
            "left_flanking_gene": self.left_flanking_gene,
            "right_flanking_gene": self.right_flanking_gene,
            "segment_color": self.segment_color
        }


    def add_node(self, seg, index, color, x, y, sample):
        node_id = seg.split(" (")[0]
        td = self.tooltip_data.get((sample, node_id), {})
        self.nodes_dict[index] = {
            "id": seg,
            "index": index,
            "color": color,
            "radius": self.node_radius,
            "full_name": td.get("full_name", node_id),
            "strand": td.get("strand", ""),
            "status": td.get("status", ""),
            "SNPs": td.get("SNPs", ""),
            "Insertions": td.get("Insertions", ""),
            "Deletions": td.get("Deletions", ""),
            "function_segment": td.get("function_segment", ""),
            "found_rss_count": td.get("found_rss_count", ""),
            "5_RSS_seq": td.get("5_RSS_seq", ""),
            "3_RSS_seq": td.get("3_RSS_seq", ""),
            "x": x,
            "y": y
        }
        self.node_ids.add(node_id)


    def add_common_node(self, index, seg, x, hap1_prev, hap2_prev):
        node_id = seg.split(" (")[0]
        sample = self.accessions[0]

        color = self.segment_color.get(self.tooltip_data.get((sample, node_id), {}).get("segment", "N/A"), "white")
        index += 1
        self.add_node(seg, index, color, x, self.y_positions["common"], sample)
        self.add_link(hap1_prev, index)
        self.add_link(hap2_prev, index)
        return index


    def add_haplotype_node(self, index, seg, x, haplotype, hap_prev):
        node_id = seg.split(" (")[0]
        if haplotype == "haplotype1":
            sample = self.accessions[0]
        else:
            sample = self.accessions[1]

        color = self.segment_color.get(self.tooltip_data.get((sample, node_id), {}).get("segment", "N/A"), "white")
        index += 1
        self.add_node(seg, index, color, x, self.y_positions[haplotype], sample)
        self.add_link(hap_prev, index)
        return index


    def add_link(self, source, target):
        self.links.append({"source": str(source), "target": str(target)})
