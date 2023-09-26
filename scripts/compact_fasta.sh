#!/bin/bash
: <<'END_COMMENT'
Bash script for condensing the default fasta format to one line.
For all fasta files in TR/DB and TR/online.
END_COMMENT

current_dir=$(cd .. && pwd)"

if [ -d "${current_dir}/TR" ]; then
    for dir in "${current_dir}/TR/"*; do
        if [ -d "${dir}" ]; then
            for fasta_file in "${dir}"/*; do
                if [ -f "${fasta_file}" ]; then
                    compact_dir="${dir}/compact"
                    if [ ! -d "${compact_dir}" ]; then
                        mkdir "${compact_dir}"
                    fi
                    output_file="${compact_dir}/compact_$(basename ${fasta_file})"
                    awk '/^>/ {if (NR>1) printf("\n"); printf("%s@", $0); next} {printf("%s", $0)} END {printf("\n")}' "${fasta_file}" > "${output_file}"
                fi
            done
        fi
    done
fi 
