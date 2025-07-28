#!/usr/bin/env python3
import os
import re
import pandas as pd
import subprocess
from Bio import SeqIO, AlignIO
from Bio.Seq import Seq
import matplotlib.pyplot as plt
import logomaker

# ─── CONFIGURATION ─────────────────────────────────────────────────────────
#BASE_DIR = "/Users/ott/PycharmProjects/VDJ-insights/Broken_regions/rhesus/"
BASE_DIR = "/path/to/outcomes/HPGC/"
# Annotation files
bcr = pd.read_excel(os.path.join(BASE_DIR, "scaffolds", "bcr", "annotation", "annotation_report_all.xlsx"))
tcr = pd.read_excel(os.path.join(BASE_DIR, "scaffolds", "tcr", "annotation", "annotation_report_all.xlsx"))

# FASTA inputs
FASTA_B = os.path.join(BASE_DIR, "scaffolds", "bcr", "tmp", "CDR", "blast", "report", "target_sequence.fasta")
FASTA_T = os.path.join(BASE_DIR, "scaffolds", "tcr", "tmp", "CDR", "blast", "report", "target_sequence.fasta")

# Directory for CDR MSAs and output
ALIGN_DIR = os.path.join(BASE_DIR, "figure", "cdr_logos", "cdr_alignments")
LOGO_DIR = os.path.join(BASE_DIR, "figure", "cdr_logos")
OUTPUT_XLSX = os.path.join(LOGO_DIR, "cdr_extracted.xlsx")
CLUSTALO = "/path/to/miniconda3/envs/bin/clustalo"

os.makedirs(LOGO_DIR, exist_ok=True)
os.makedirs(ALIGN_DIR, exist_ok=True)

# ─── GLOBAL COLOR MAPPING ───────────────────────────────────────────────────
# Use a fixed, colorblind-friendly palette for all amino acids
AMINO_ACIDS = [
    'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L',
    'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y', '*'
]
#AMINO_ACIDS = ['A', 'T', 'G', 'C']   bps logos

cmap = plt.get_cmap('tab20')
colors = cmap(range(len(AMINO_ACIDS)))
# Convert RGBA tuples to hex strings for logomaker compatibility
COLOR_DICT = {aa: colors[k] for k, aa in enumerate(AMINO_ACIDS)}

# Grey for gaps
COLOR_DICT['-'] = '#CCCCCC'

# ─── HELPERS ────────────────────────────────────────────────────────────────
def parse_short_name(name):
    m = re.match(r"(IG[HKL]|TR[ABDG])([VDJ])", name)
    return (m.group(1), m.group(2)) if m else (None, None)

def run(cmd):
    print(f"→ Running: {' '.join(cmd)}")
    subprocess.run(cmd, check=True)

def extract_cdrs(df, seg_fasta):
    fasta_dict = {rec.id: rec.seq for rec in SeqIO.parse(seg_fasta, "fasta")}
    df = df.dropna(subset=["CDR1_start", "CDR1_stop", "CDR2_start", "CDR2_stop"]).copy()
    df[["Locus", "Class"]] = df["Short name"].apply(lambda x: pd.Series(parse_short_name(x)))
    records = []

    for idx, row in df.iterrows():
        short = row["Short name"]
        sample = row.get("Sample", "")
        region = row["Locus"]
        clas = row["Class"]
        key = f"{short}_{idx}"
        seq = fasta_dict.get(key)
        if seq is None:
            print(f"⚠️  No FASTA record for {key}, skipping.")
            continue
        if row.get("Strand") == "-":
            seq = seq.reverse_complement()

        c1_start = int(row["CDR1_start"])
        c1_stop = int(row["CDR1_stop"])
        c2_start = int(row["CDR2_start"])
        c2_stop = int(row["CDR2_stop"])

        cdr1_nt = str(seq[c1_start:c1_stop])
        cdr2_nt = str(seq[c2_start:c2_stop])
        cdr1_aa = str(Seq(cdr1_nt).translate(to_stop=False))
        cdr2_aa = str(Seq(cdr2_nt).translate(to_stop=False))

        records.append({
            "Sample": sample,
            "Region": region,
            "Segment_type": clas,
            "Short name": short,
            "Index": idx,
            "FASTA key": key,
            "CDR1_start": c1_start,
            "CDR1_stop": c1_stop,
            "CDR1_bps": cdr1_nt,
            "CDR1_aa": cdr1_aa,
            "CDR2_start": c2_start,
            "CDR2_stop": c2_stop,
            "CDR2_bps": cdr2_nt,
            "CDR2_aa": cdr2_aa,
        })

    return pd.DataFrame(records)

# ─── MAIN WORKFLOW ─────────────────────────────────────────────────────────
def main():
    # Load or extract CDR table
    if os.path.exists(OUTPUT_XLSX):
        out_df = pd.read_excel(OUTPUT_XLSX)
    else:
        df_b = extract_cdrs(bcr, FASTA_B)
        df_t = extract_cdrs(tcr, FASTA_T)
        out_df = pd.concat([df_b, df_t], ignore_index=True)
        out_df.to_excel(OUTPUT_XLSX, index=False)
    out_df.drop_duplicates('Short name', inplace=True)

    # Compute number of sequences per region
    region_counts = out_df.groupby('Region').size().to_dict()

    # Gather PWM data
    logo_data = {}
    for locus, grp in out_df.groupby('Region'):
        logo_data[locus] = {}
        for cdr in ('CDR1','CDR2'):
            seqs = grp[f'{cdr}_aa'].dropna().tolist()
            #seqs = grp[f'{cdr}_bps'].dropna().tolist() #bps logos
            if not seqs:
                continue
            fasta_f = os.path.join(ALIGN_DIR, f"{locus}_{cdr}.fasta")
            aln_f = os.path.join(ALIGN_DIR, f"{locus}_{cdr}.aln")
            # Write & align
            with open(fasta_f, 'w') as fh:
                for i, aa in enumerate(seqs):
                    fh.write(f">{locus}_{cdr}_{i}\n{aa}\n")
            run([CLUSTALO, '-i', fasta_f, '-o', aln_f, '--force', '--outfmt', 'clustal'])
            aln = AlignIO.read(aln_f, 'clustal')
            df_aln = pd.DataFrame([list(rec.seq) for rec in aln])
            freqs = {col: df_aln[col].value_counts(normalize=True).to_dict()
                      for col in df_aln.columns}
            pwm = pd.DataFrame(freqs).fillna(0).T
            if pwm.empty:
                continue
            # Separate gap frequencies
            gap = pwm.get('-', pd.Series(0, index=pwm.index))
            aa_pwm = pwm.drop('-', axis=1, errors='ignore')
            logo_data[locus][cdr] = (aa_pwm, gap)

    # Plot combined grid with colorblind-friendly palette + counts
    regions = sorted(logo_data.keys())
    n = len(regions)
    fig, axs = plt.subplots(n,2,figsize=(12, 3 * n),squeeze=False)
    # Column titles
    axs[0, 0].set_title('CDR1',fontsize=14)
    axs[0, 1].set_title('CDR2',fontsize=14)

    for i, locus in enumerate(regions):
        count = region_counts.get(locus,0)
        for j, cdr in enumerate(('CDR1', 'CDR2')):
            ax = axs[i, j]
            data = logo_data[locus].get(cdr)
            if not data:
                ax.axis('off'); continue
            pwm_aa, gap = data
            filtered_colors = {aa: COLOR_DICT[aa] for aa in pwm_aa.columns if aa in COLOR_DICT}
            filtered_colors['-'] = '#CCCCCC'
            logomaker.Logo(pwm_aa, ax=ax, stack_order='small_on_top', color_scheme=filtered_colors)
            trans = ax.transData + ax.get_yaxis_transform()
            for k, val in enumerate(gap): ax.vlines(k, 1, 1 + val, transform=trans, colors='gray', lw=2)
            ax.set_xticks([]); ax.set_yticks([])
            if j == 0:
                ax.set_ylabel(f"{locus}\n(N={count})", rotation=0, labelpad=30, va='center', fontsize=14)

    fig.tight_layout(rect=[0, 0, 1, 0.96])
    #fig.suptitle('CDR sequence logos per region', x=0.5, y=0.99, fontsize=16)
    out_file = os.path.join(LOGO_DIR, 'combined_cdr_logos_with_counts.svg')
    fig.savefig(out_file, format='svg')
    plt.close(fig)
    print(f"✅ Saved combined figure: {out_file}")

if __name__ == '__main__':
    main()