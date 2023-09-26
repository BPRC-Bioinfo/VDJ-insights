#!/bin/bash
: <<'END_COMMENT'
Bash Script for counting the amount of V, D or J gene sequences in the fasta file from TR/DB and TR/online.
END_COMMENT

current_dir="$(pwd)/data"
if [ -d "${current_dir}/TR/" ]; then
    for dir in "${current_dir}/TR/"*; do
        echo "${dir}" 
        if [ -d "${dir}" ]; then
            for fasta_file in "${dir}/"*; do
                if [ -f "${fasta_file}" ]; then
                    seq=$(cat "${fasta_file}" | egrep "^>" | wc -l)
                    echo "Amount of sequences: ${seq} in $(basename ${fasta_file})"
                fi
            done
        fi
    done
fi