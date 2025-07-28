#This script runs Digger for X number of samples, for IG, TR or Both regions.
#The outcomes of this script are csv files containing Digger outcomes from the regions in forward direction, reverse direction and both files combined.
#A log file is created for the Digger

import subprocess
import os
import time
import glob
import pandas as pd
import yaml
import shutil

#If you have not done so, install Digger and seqkit + read the Before_starting_Digger.txt file.

# Load the configuration file. This file includes all directories
with open("config_human.yaml", "r") as file:
    config = yaml.safe_load(file)

# Assign variables from the config
species = config["species"] #species used for motif files + how you want to name your files
Latin_species = config["Latin_species"] #species used to extract IMGT gapped segments
Latin_species_ = Latin_species.replace(" ", "_") #For the IMGT commands a space is needed, but in the created files this is a "_"
samples = config["samples"] #The list of animals included. This should be a part of the contig file names.
IG_library = config["IG_library"]
TR_library = config["TR_library"]
output_dir = config["output_dir"] #Output dir for libraries
base_region_dir = config["base_region_dir"] #Directory including all contigs which need to be annotated
digger_output_dir = config["digger_output_dir"] #Digger output folder
region = config["region"] #IG, TR or Both

### Start code!

# Ensure output directory exists
os.makedirs(output_dir, exist_ok=True)
os.makedirs(digger_output_dir, exist_ok=True)

# List of loci to extract form the library
loci_list_IG = ["IGHV", "IGHD", "IGHJ", "IGKV", "IGKJ", "IGLV", "IGLJ"]
loci_list_TR = ["TRAV", "TRAJ","TRBV","TRBD","TRBJ","TRGV", "TRGJ","TRDV", "TRDD", "TRDJ"]
loci_list_both = ["IGHV", "IGHD", "IGHJ", "IGKV", "IGKJ", "IGLV", "IGLJ","TRAV", "TRAJ","TRBV","TRBD","TRBJ","TRGV", "TRGJ","TRDV", "TRDD", "TRDJ"]

if region == "IG":
    loci_list = loci_list_IG
elif region == "TR":
    loci_list = loci_list_TR
else:
    loci_list = loci_list_both

# Extract sequences libraries using seqkit
for locus in loci_list:
    output_fasta = f"{output_dir}{species}_lib_{locus}.fasta"

    if os.path.exists(output_fasta) and os.path.getsize(output_fasta) > 0:
        print(f"✅ {output_fasta} already exists. Skipping extraction.")
        continue

    if region == "Both":
        if locus[0:2] == "IG":
            command = f"seqkit grep -r -p '{locus}' {IG_library} >> {output_fasta}"
        else:
            command = f"seqkit grep -r -p '{locus}' {TR_library} >> {output_fasta}"
    elif region == "IG":
        command = f"seqkit grep -r -p '{locus}' {IG_library} -o {output_fasta}"
    else:  # region == "TR"
        command = f"seqkit grep -r -p '{locus}' {TR_library} -o {output_fasta}"
    subprocess.run(command, shell=True, check=True)
    print(f"Extracted {locus} sequences into {output_fasta}")

'''
#Combine libraries, is needed for digger -ref imgt,{species}_lib_{region}V(D)J.fasta
#For now we have decided not to include this part in the code, since it gives errors of IGK, IGL and maybe more...
#Extract if region == "IG" or "both":
list_cat =[]

if region == "IG" or region == "Both":
    cat_IGH = f"{output_dir}{species}_lib_IGHVDJ.fasta"
    cat_IGK = f"{output_dir}{species}_lib_IGKVJ.fasta"
    cat_IGL = f"{output_dir}{species}_lib_IGLVJ.fasta"
    list_cat += [cat_IGH, cat_IGK, cat_IGL]  # Make sure list is created properly

if region == "TR" or region == "Both":
    cat_TRA = f"{output_dir}{species}_lib_TRAVJ.fasta"
    cat_TRB = f"{output_dir}{species}_lib_TRBVDJ.fasta"
    cat_TRD = f"{output_dir}{species}_lib_TRDVDJ.fasta"
    cat_TRG = f"{output_dir}{species}_lib_TRGVJ.fasta"
    list_cat += [cat_TRA, cat_TRB, cat_TRD, cat_TRG]  # Make sure list is created properly
'''
for locus in loci_list:
    '''
    for i in range(len(list_cat)):
        current_list = list_cat[i]
        if os.path.exists(current_list) and os.path.getsize(current_list) > 0:
            continue
        else:
           # Expand the wildcard and get all matching files
           locus_ = locus[0:3]
           command = f"cat {output_dir}{species}_lib_{locus_}*.fasta > {current_list}"
           result = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True, text=True)
           if result.returncode == 0:
               print(f"✅ Successfully combined libraries into {current_list}")
    '''
#Extract gapped libraries from the IMGT website
    gap_lib_end = f"{output_dir}{Latin_species_}_{locus[0:3]}V_gapped.fasta"
    if os.path.exists(gap_lib_end) and os.path.getsize(gap_lib_end) > 0:
        continue

    gap_lib = f"{Latin_species_}_{locus[0:3]}V_gapped.fasta"
    gap_fix = f"{Latin_species_}_{locus[0:3]}V_gapped_fixed.fasta"
    extract_command = ['extract_refs', '-L', locus[0:3], f'"{Latin_species}"']
    subprocess.run(' '.join(extract_command),shell=True, check=True)
    if species == "rhesus_macaque" and locus[0:2] == "IG": #Digger fixes gaps for the rhesus macaque, but does not do this for humans...
        fix_gaps_command = ["fix_macaque_gaps", gap_lib, gap_fix, locus[0:3]]
        subprocess.run(' '.join(fix_gaps_command),shell=True, check=True)
        shutil.copyfile(gap_fix, gap_lib_end)
    else:
        shutil.copyfile(gap_lib,gap_lib_end)


# Gives dictionary of all files per region, will be the input for the digger run
loci_info = {
    "IGH": {
        "v_ref": f"{output_dir}{species}_lib_IGHV.fasta",
        "d_ref": f"{output_dir}{species}_lib_IGHD.fasta",
        "j_ref": f"{output_dir}{species}_lib_IGHJ.fasta",
        "v_ref_gap": f"{output_dir}{Latin_species_}_IGHV_gapped.fasta",
        "combined": f"{output_dir}{species}_lib_IGHVDJ.fasta"},
    "IGK": {
        "v_ref": f"{output_dir}{species}_lib_IGKV.fasta",
        "j_ref": f"{output_dir}{species}_lib_IGKJ.fasta",
        "v_ref_gap": f"{output_dir}{Latin_species_}_IGKV_gapped.fasta",
        "combined": f"{output_dir}{species}_lib_IGKVJ.fasta"},
    "IGL": {
        "v_ref": f"{output_dir}{species}_lib_IGLV.fasta",
        "j_ref": f"{output_dir}{species}_lib_IGLJ.fasta",
        "v_ref_gap": f"{output_dir}{Latin_species_}_IGLV_gapped.fasta",
        "combined": f"{output_dir}{species}_lib_IGLVJ.fasta"},
    "TRA": {
        "v_ref": f"{output_dir}{species}_lib_TRAV.fasta",
        "j_ref": f"{output_dir}{species}_lib_TRAJ.fasta",
        "v_ref_gap": f"{output_dir}{Latin_species_}_TRAV_gapped.fasta",
        "combined": f"{output_dir}{species}_lib_TRAVJ.fasta"},
    "TRB": {
        "v_ref": f"{output_dir}{species}_lib_TRBV.fasta",
        "d_ref": f"{output_dir}{species}_lib_TRBD.fasta",
        "j_ref": f"{output_dir}{species}_lib_TRBJ.fasta",
        "v_ref_gap": f"{output_dir}{Latin_species_}_TRBV_gapped.fasta",
        "combined": f"{output_dir}{species}_lib_TRBVDJ.fasta"},
    "TRG": {
        "v_ref": f"{output_dir}{species}_lib_TRGV.fasta",
        "j_ref": f"{output_dir}{species}_lib_TRGJ.fasta",
        "v_ref_gap": f"{output_dir}{Latin_species_}_TRGV_gapped.fasta",
        "combined": f"{output_dir}{species}_lib_TRGVJ.fasta"},
    "TRD": {
        "v_ref": f"{output_dir}{species}_lib_TRDV.fasta",
        "d_ref": f"{output_dir}{species}_lib_TRDD.fasta",
        "j_ref": f"{output_dir}{species}_lib_TRDJ.fasta",
        "v_ref_gap": f"{output_dir}{Latin_species_}_TRDV_gapped.fasta",
        "combined": f"{output_dir}{species}_lib_TRDVDJ.fasta"}
}

# Define digger log file
log_file = f"{digger_output_dir}digger_run.log"

def run_command(command, output_csv, log_file):
    """Runs a command and logs output/errors."""
    with open(log_file, "a") as log:
        log.write(f"\nRunning: {' '.join(command)}\n")  # Log command
        process = subprocess.run(command, stdout=log, stderr=log, text=True)

    return process


# Run Digger for each locus
for sample in samples:
    for locus, refs in loci_info.items():
        if locus in ["TRA", "TRD"]: #In case of TRA and TRD the same region is taken.
            region_fasta_files = glob.glob(f"{base_region_dir}/{sample}*_TRA.fasta")
        else:
            region_fasta_files = glob.glob(f"{base_region_dir}/{sample}*_{locus}.fasta")

        if not region_fasta_files: #Result in no runs for regions which are not selected
            if (region == "IG" and locus.startswith("TR")) or (region == "TR" and locus.startswith("IG")):
                continue
            print(f"⚠️ No region file found for {locus} in {sample}. Skipping Digger run.")
            continue

        region_fasta = region_fasta_files[0]
        output_csv_forwards = f"{digger_output_dir}{sample}_{locus}_digger_forwards.csv"
        output_csv_reverse = f"{digger_output_dir}{sample}_{locus}_digger_reverse.csv"
        combined_csv = f"{digger_output_dir}{sample}_{locus}_digger_combined.csv"

        # Skip if both output files exist
        if os.path.exists(output_csv_forwards) and os.path.exists(output_csv_reverse) and os.path.exists(combined_csv):
            print(f"✅ Digger output for {locus} in {sample} already exists. Skipping.")
            continue

        # Construct Digger commands, both forward and reverse
        if locus in ["IGH", "TRB", "TRD"]:
            command_f = ["digger", region_fasta, "-locus", locus, "-species", species,
                         "-v_ref", refs["v_ref"], "-d_ref", refs["d_ref"], "-j_ref", refs["j_ref"],
                         "-v_ref_gapped", refs["v_ref_gap"],"-sense", "forward",
                         output_csv_forwards]
            command_r = command_f[:-3] + ["-sense", "reverse", output_csv_reverse]
        else:
            command_f = ["digger", region_fasta, "-locus", locus, "-species", species,
                         "-v_ref", refs["v_ref"], "-j_ref", refs["j_ref"],
                         "-v_ref_gapped", refs["v_ref_gap"], "-sense", "forward",
                         output_csv_forwards]
            command_r = command_f[:-3] + ["-sense", "reverse", output_csv_reverse]

        # Run forward command and log output
        start_time_f = time.perf_counter()
        process_f = run_command(command_f, output_csv_forwards, log_file)
        end_time_f = time.perf_counter()
        time.sleep(2) #Digger needs some time to shut down, also to prevent digger to delete files which are still in use in the previous run.

        if process_f.returncode == 0:
            print(f"✅ Digger ran {locus} of {sample} (forwards) in {end_time_f - start_time_f:.2f}s.")
        else:
            print(f"❌ Error running Digger (forwards) for {locus} of {sample}. Check log: {log_file}")

        # Run reverse command and log output
        start_time_r = time.perf_counter()
        process_r = run_command(command_r, output_csv_reverse, log_file)
        end_time_r = time.perf_counter()
        time.sleep(2)

        if process_r.returncode == 0:
            print(f"✅ Digger ran {locus} of {sample} (reverse) in {end_time_r - start_time_r:.2f}s.")
        else:
            print(f"❌ Error running Digger (reverse) for {locus} of {sample}. Check log: {log_file}")

        #This part combines the forward and reverse csv outcome to a combined file.
        try:
            df_f = pd.read_csv(output_csv_forwards)
            df_r = pd.read_csv(output_csv_reverse)
            combined_df = pd.concat([df_f, df_r], ignore_index=True)
            combined_df.to_csv(combined_csv, index=False)
            print(f"✅ Successfully merged forward and reverse Digger results into {combined_csv}")

        except Exception as e:
            print(f"❌ Error merging CSV files for {locus} in {sample}: {e}")

#Digger creates many files which are not usefull after the run. Here I try to clean already a little bit.
for file in glob.glob(f"*{Latin_species_}_*"):  # Matches all files starting with "*Gorilla_gorilla*"
    os.remove(file)

for file in glob.glob(f"*{species}_*"):  # Matches all files starting with "*gorilla_*"
    os.remove(file)

print(f"✅ Everything is cleaned up!")

