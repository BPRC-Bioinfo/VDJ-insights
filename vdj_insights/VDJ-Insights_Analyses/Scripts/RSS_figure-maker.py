import os
import re
import pandas as pd
import matplotlib.pyplot as plt
import logomaker
from Bio import SeqIO
from matplotlib.lines import Line2D

# ─── Helpers ───────────────────────────────────────────────────────────────

def parse_short_name(name):
    m = re.match(r"(IG[HKL]|TR[ABDG])([VDJ])", name)
    return (m.group(1), m.group(2)) if m else (None, None)

def split_rss(seq, side):
    s = str(seq)
    return (s[:7], s[-9:]) if side=="hept-first" else (s[:9], s[-7:])
    
def make_matrix(seqs, length):
    filtered = [s for s in seqs if s and len(s)==length]
    return logomaker.alignment_to_matrix(filtered, to_type="counts") if filtered else None

# ─── Load & Parse ─────────────────────────────────────────────────────────

base = "/path/to/outcomes/HPGC/"
bcr = pd.read_excel(os.path.join(base, "scaffolds","bcr","annotation","annotation_report_all.xlsx"))
tcr = pd.read_excel(os.path.join(base, "scaffolds","tcr","annotation","annotation_report_all.xlsx"))
df  = pd.concat([bcr, tcr], ignore_index=True)

df[["Locus","Class"]] = df["Short name"].apply(lambda x: pd.Series(parse_short_name(x)))

# ─── Extract all RSS sequences from FASTA ──────────────────────────────────
d_rows = df[df["Class"] == "D"].copy()
non_d  = df[df["Class"] != "D"].copy()

d5 = d_rows.copy(); d5["RSS_side"] = "5"
d3 = d_rows.copy(); d3["RSS_side"] = "3"
df = pd.concat([non_d, d5, d3], ignore_index=True)

# 2) V and J get their side by class
df.loc[df["Class"] == "V", "RSS_side"] = "3"   # V always 3′-flank
df.loc[df["Class"] == "J", "RSS_side"] = "5"   # J always 5′-flank

_fasta_cache = {}
def get_seq_record(path):
    if path not in _fasta_cache:
        _fasta_cache[path] = SeqIO.read(path, "fasta")
    return _fasta_cache[path]
    
def extract_rss_by_coords(row):
    rec    = get_seq_record(row["Path"])
    Locus  = row["Locus"]
    seg    = row["Class"]
    side   = row["RSS_side"]
    strand = row["Strand"]
    start  = int(row["Start coord"])
    end    = int(row["End coord"])
    seqlen = len(rec.seq)

    # determine spacer length and region total length
    if seg == "V":
        spacer    = 23 if Locus in ("IGH","IGL","TRA","TRB","TRD","TRG") else 12
    elif seg == "J":
        spacer    = 12 if Locus in ("IGL","TRA","TRB","TRD","TRG") else 23
    else:  # D flank
        if side == "5":
            spacer = 12
        else:  # '3'
            spacer = 12 if Locus == "IGH" else 23

    region_len = 7 + spacer + 9  # heptamer + spacer + nonamer

    # choose genomic window around gene boundary
    if side == "3":
        # flank on the 3′ end of the gene
        if strand == "+":
            win_start = end
            win_end   = end + region_len
        else:
            win_end   = start
            win_start = start - region_len
    else:  # side == "5"
        # flank on the 5′ end of the gene
        if strand == "+":
            win_end   = start
            win_start = start - region_len
        else:
            win_start = end
            win_end   = end + region_len

    # clamp to scaffold boundaries
    win_start = max(0, win_start)
    win_end   = min(seqlen, win_end)

    raw = rec.seq[win_start:win_end]
    # orient to gene (+) so that heptamer always at left for side="3", at right for side="5"
    oriented = raw if strand=="+" else raw.reverse_complement()
    return str(oriented)
            
df["RSS_seq"] = df.apply(extract_rss_by_coords, axis=1)

idx = ["Sample", "Region", "Locus", "Class", "Start coord", "End coord", "Short name", "Function", "Strand"]

df = (df.pivot_table(index=idx, columns="RSS_side", values="RSS_seq", aggfunc="first").reset_index())
df = df.rename(columns={"5": "RSS_seq_5", "3": "RSS_seq_3"})

out_path = os.path.join(base, "figure", "rss_regions_all.xlsx")
export_df = df[["Sample", "Region", "Locus", "Class", "Start coord", "End coord", "Strand", "Short name", "Function", "RSS_seq_5", "RSS_seq_3"]]
export_df.to_excel(out_path, index=False)
print(f"Wrote RSS list to {out_path}")
print("Number of RSS:", len(df))

# ─── Build logo_data ──────────────────────────────────────────────────────
#Drops Segments mapping on wrong locus (especially done for D segments)
mask = ((df["Locus"] == df["Region"]) | ((df["Region"] == "TRA") & (df["Locus"] == "TRD")))
df = df[mask].reset_index(drop=True)
print("Number of RSS, minus wrong-region:", len(df))
#Drop RSS with NNNs
df = df.loc[~df["RSS_seq_5"].str.contains("N", na=False) & ~df["RSS_seq_3"].str.contains("N", na=False)].reset_index(drop=True)
print("Number of RSS, minus wrong-region and NNN:", len(df))

df.drop_duplicates(subset=["Short name", "RSS_seq_5", "RSS_seq_5"], keep='first', inplace=True)
print("Number of RSS, minus duplicates, wrong-region and NNN:", len(df))

counts = df.groupby(['Locus', 'Class']).size()
region_counts = {f"{locus}{cls}": int(cnt) for (locus, cls), cnt in counts.items()}
df = df.dropna(subset=["Locus","Class","Strand"])

logo_data = {}
for locus in sorted(df["Locus"].unique()):
    V_h, V_n = [], []
    J_h, J_n = [], []
    D5_h, D5_n = [], []
    D3_h, D3_n = [], []
    
    for _, row in df[df["Locus"]==locus].iterrows():
        seg = row["Class"]

        if seg=="V":
            seq = row["RSS_seq_3"]
            h,n = split_rss(seq, "hept-first")
            V_h.append(h); V_n.append(n)
        elif seg=="J":
            seq = row["RSS_seq_5"]
            n,h = split_rss(seq, "non-first")
            J_n.append(n); J_h.append(h)
        elif seg=="D":
            seq5 = row["RSS_seq_5"]
            seq3 = row["RSS_seq_3"]
            if pd.notna(seq5):
                n5, h5 = split_rss(seq5, "non-first")
                D5_n.append(n5); D5_h.append(h5)
            if pd.notna(seq3):
                h3, n3 = split_rss(seq3, "hept-first")
                D3_h.append(h3); D3_n.append(n3)

    mats = {
        "V_hept": make_matrix(V_h, 7),
        "V_non":  make_matrix(V_n, 9),
        "J_non":  make_matrix(J_n, 9),
        "J_hept": make_matrix(J_h, 7),
    }
    if locus in ("IGH","TRB","TRD"):
        mats.update({
            "D5_non":  make_matrix(D5_n, 9),
            "D5_hept": make_matrix(D5_h, 7),
            "D3_hept": make_matrix(D3_h, 7),
            "D3_non":  make_matrix(D3_n, 9),
        })


    logo_data[locus] = mats

# ─── Plotting ─────────────────────────────────────────────────────────────

vj_loci  = ["IGH","IGK","IGL","TRA","TRB","TRD","TRG"]
d_loci   = ["IGH","TRB","TRD"]
row_specs= [(l,"VJ") for l in vj_loci] + [(l,"D") for l in d_loci]

vj_parts   = ["V_hept","V_non","J_non","J_hept"]
vj_labels  = ["V-Heptamer","V-Nonamer","J-Nonamer","J-Heptamer"]
d_parts    = ["D5_non","D5_hept","D3_hept","D3_non"]
d_labels   = ["5′D-Nonamer","5′D-Heptamer","3′D-Heptamer","3′D-Nonamer"]

n_rows = len(row_specs)
fig, axs = plt.subplots(
    n_rows, 4,
    figsize=(12, 2.5*n_rows),
    gridspec_kw={"wspace":0.5,"hspace":0.6}
)

# Column titles
for j,lbl in enumerate(vj_labels):
    axs[0,j].set_title(lbl, fontsize=14, pad=12)
for j,lbl in enumerate(d_labels):
    axs[len(vj_loci), j].set_title(lbl, fontsize=14, pad=12)

# Draw each row
for i,(locus,segtype) in enumerate(row_specs):
    count_v = region_counts.get(locus + "V", 0)
    count_j = region_counts.get(locus + "J", 0)
    parts = vj_parts if segtype=="VJ" else d_parts

    # 1) Logos
    for j,part in enumerate(parts):
        ax=axs[i,j]
        mat=logo_data[locus].get(part)
        if mat is None: ax.axis("off")
        else: logomaker.Logo(mat, ax=ax, show_spines=False, color_scheme='classic')
        ax.set_xticks([]); ax.set_yticks([])

    # 2) Thin spacer‐bars
    p0,p1,p2,p3 = [axs[i,c].get_position() for c in (0,1,2,3)]
    y_mid = (p0.y0 + p0.y1)/2

    def draw_bar(x0,x1,text):
        bar = Line2D([x0,x1],[y_mid,y_mid], transform=fig.transFigure, lw=1.2)
        fig.add_artist(bar)
        fig.text((x0+x1)/2, y_mid+0.01, text, ha="center", va="bottom", fontsize=9)

    # determine spacers exactly as you specified
    if segtype=="VJ":
        v_bp = 23 if locus in ["IGH","IGL","TRA","TRB","TRD","TRG"] else 12
        j_bp = 12 if locus in ["IGL","TRA","TRB","TRD","TRG"] else 23
    else:  # D
        v_bp = 12
        j_bp = 12 if locus=="IGH" else 23

    # V‐side
    raw0,raw1 = p0.x1, p1.x0
    inset = (raw1-raw0)*0.2
    draw_bar(raw0+inset, raw1-inset, f"{v_bp} bp")

    # J‐side (or downstream D)
    raw2,raw3 = p2.x1, p3.x0
    inset2 = (raw3-raw2)*0.2
    draw_bar(raw2+inset2, raw3-inset2, f"{j_bp} bp")

    # 3) Thick D‐connector + D‐label
    if segtype=="D":
        count_d = region_counts.get(locus + "D", 0)
        raw_l,raw_r = p1.x1, p2.x0
        gap = raw_r-raw_l; di=gap*0.1
        xs,xe = raw_l+di, raw_r-di
        fig.add_artist(Line2D([xs,xe],[y_mid,y_mid], transform=fig.transFigure, lw=3.0))
        fig.text((xs+xe)/2, y_mid+0.01, f"{locus}D\n(N = {count_d})",
                 ha="center", va="bottom", fontsize=12)

    # 4) Dual row‐labels
    left_lbl  = locus + ("V" if segtype=="VJ" else "D") + f"\n(N={count_v})"
    right_lbl = locus + ("J" if segtype=="VJ" else "D") + f"\n(N={count_j})"
    axL,axR = axs[i,0], axs[i,-1]
    axL.yaxis.set_label_position("left");  axL.yaxis.tick_left()
    axR.yaxis.set_label_position("right"); axR.yaxis.tick_right()
    if segtype!="D":
        axL.set_ylabel(left_lbl,  rotation=0, labelpad=30, va="center", ha="right", fontsize=12)
        axR.set_ylabel(right_lbl, rotation=0, labelpad=30, va="center", ha="left",  fontsize=12)

# Title & save
fig.subplots_adjust(top=0.88)
fig.suptitle("RSS sequence logos across the IG & TCR loci", y=0.92, fontsize=18)

out = os.path.join(base, "figure", "rss-figure")
fig.savefig(f"{out}.svg", format="svg", bbox_inches="tight")
fig.savefig(f"{out}.pdf", format="pdf", bbox_inches="tight")
plt.close(fig)

# ─── Second figure: only Functional segments ───────────────────────────────

# 1) filter
df_func = df[df["Function"] == "Functional"].reset_index(drop=True)

# 2) rebuild logo_data for this subset
counts_f = df_func.groupby(["Locus","Class"]).size()
region_counts_f = {
    f"{locus}{cls}": int(cnt)
    for (locus, cls), cnt in counts_f.items()
}

logo_data_f = {}
for locus in sorted(df_func["Locus"].unique()):
    V_h, V_n = [], []
    J_h, J_n = [], []
    D5_h, D5_n = [], []
    D3_h, D3_n = [], []

    sub = df_func[df_func["Locus"] == locus]
    for _, row in sub.iterrows():
        seg = row["Class"]
        if seg == "V":
            seq = row["RSS_seq_3"]
            h,n = split_rss(seq, "hept-first")
            V_h.append(h); V_n.append(n)
        elif seg == "J":
            seq = row["RSS_seq_5"]
            n,h = split_rss(seq, "non-first")
            J_n.append(n); J_h.append(h)
        elif seg == "D":
            if pd.notna(row["RSS_seq_5"]):
                n5,h5 = split_rss(row["RSS_seq_5"], "non-first")
                D5_n.append(n5); D5_h.append(h5)
            if pd.notna(row["RSS_seq_3"]):
                h3,n3 = split_rss(row["RSS_seq_3"], "hept-first")
                D3_h.append(h3); D3_n.append(n3)

    def make_matrix(seqs, length):
        good = [s for s in seqs if len(s)==length]
        return logomaker.alignment_to_matrix(good, to_type="counts") if good else None

    mats = {
        "V_hept": make_matrix(V_h, 7),
        "V_non":  make_matrix(V_n, 9),
        "J_non":  make_matrix(J_n, 9),
        "J_hept": make_matrix(J_h, 7),
    }
    if locus in ("IGH","TRB","TRD"):
        mats.update({
            "D5_non":  make_matrix(D5_n, 9),
            "D5_hept": make_matrix(D5_h, 7),
            "D3_hept": make_matrix(D3_h, 7),
            "D3_non":  make_matrix(D3_n, 9),
        })
    logo_data_f[locus] = mats

# 3) plotting (exactly same as before, just substitute logo_data_f & region_counts_f)
fig, axs = plt.subplots(
    n_rows, 4,
    figsize=(12, 2.5*n_rows),
    gridspec_kw={"wspace":0.5,"hspace":0.6}
)

# Column titles
for j,lbl in enumerate(vj_labels):
    axs[0,j].set_title(lbl, fontsize=14, pad=12)
for j,lbl in enumerate(d_labels):
    axs[len(vj_loci), j].set_title(lbl, fontsize=14, pad=12)

# Draw each row
for i,(locus,segtype) in enumerate(row_specs):
    count_v = region_counts_f.get(locus + "V", 0)
    count_j = region_counts_f.get(locus + "J", 0)
    parts = vj_parts if segtype=="VJ" else d_parts

    # 1) Logos
    for j,part in enumerate(parts):
        ax=axs[i,j]
        mat=logo_data_f[locus].get(part)
        if mat is None: ax.axis("off")
        else:           logomaker.Logo(mat, ax=ax, show_spines=False, color_scheme='classic')
        ax.set_xticks([]); ax.set_yticks([])

    # 2) Thin spacer‐bars
    p0,p1,p2,p3 = [axs[i,c].get_position() for c in (0,1,2,3)]
    y_mid = (p0.y0 + p0.y1)/2

    def draw_bar(x0,x1,text):
        bar = Line2D([x0,x1],[y_mid,y_mid], transform=fig.transFigure, lw=1.2)
        fig.add_artist(bar)
        fig.text((x0+x1)/2, y_mid+0.01, text, ha="center", va="bottom", fontsize=9)

    # determine spacers exactly as you specified
    if segtype=="VJ":
        v_bp = 23 if locus in ["IGH","IGL","TRA","TRB","TRD","TRG"] else 12
        j_bp = 12 if locus in ["IGL","TRA","TRB","TRD","TRG"] else 23
    else:  # D
        v_bp = 12
        j_bp = 12 if locus=="IGH" else 23

    # V‐side
    raw0,raw1 = p0.x1, p1.x0
    inset = (raw1-raw0)*0.2
    draw_bar(raw0+inset, raw1-inset, f"{v_bp} bp")

    # J‐side (or downstream D)
    raw2,raw3 = p2.x1, p3.x0
    inset2 = (raw3-raw2)*0.2
    draw_bar(raw2+inset2, raw3-inset2, f"{j_bp} bp")

    # 3) Thick D‐connector + D‐label
    if segtype=="D":
        count_d = region_counts_f.get(locus + "D", 0)
        raw_l,raw_r = p1.x1, p2.x0
        gap = raw_r-raw_l; di=gap*0.1
        xs,xe = raw_l+di, raw_r-di
        fig.add_artist(Line2D([xs,xe],[y_mid,y_mid], transform=fig.transFigure, lw=3.0))
        fig.text((xs+xe)/2, y_mid+0.01, f"{locus}D\n(N = {count_d})",
                 ha="center", va="bottom", fontsize=12)

    # 4) Dual row‐labels
    left_lbl  = locus + ("V" if segtype=="VJ" else "D") + f"\n(N={count_v})"
    right_lbl = locus + ("J" if segtype=="VJ" else "D") + f"\n(N={count_j})"
    axL,axR = axs[i,0], axs[i,-1]
    axL.yaxis.set_label_position("left");  axL.yaxis.tick_left()
    axR.yaxis.set_label_position("right"); axR.yaxis.tick_right()
    if segtype!="D":
        axL.set_ylabel(left_lbl,  rotation=0, labelpad=30, va="center", ha="right", fontsize=12)
        axR.set_ylabel(right_lbl, rotation=0, labelpad=30, va="center", ha="left",  fontsize=12)

# 4) save under a new name
out2 = os.path.join(base, "figure", "rss-figure-functional")
fig.savefig(f"{out2}.svg", format="svg", bbox_inches="tight")
fig.savefig(f"{out2}.pdf", format="pdf", bbox_inches="tight")
plt.close(fig)
print(f"Saved Functional‐only RSS logo figure to {out2}")