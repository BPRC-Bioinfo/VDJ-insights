from Bio import SeqIO
import numpy as np
import json
import gzip



def process_reads(fastq_file="data/reads/SRR8137614.fastq.gz"):
    """
    Process raw reads to get various statistics such as smallest and
    largest reads, mean read length,
    the average GC% and Q-score of the reads and the GC percentage
    and the Q-score for every position.
    """
    smallest_read = float("inf")
    largest_read = float("-inf")
    read_lens = []
    avg_q_score_read = {}
    reads_pos = {}
    reads_pos_score = {}
    gc_per_position = {}
    gc_total = 0.0
    for record in SeqIO.parse(fastq_file, "fastq"):
        len_seq = len(record.seq)
        read_lens.append(len_seq)
        largest_read = max(largest_read, len_seq)
        smallest_read = min(smallest_read, len_seq)

        q_score = record.letter_annotations["phred_quality"]
        avg_q_score_read[record.id] = sum(q_score) / len(q_score)

        for pos, score in enumerate(q_score):
            reads_pos[pos] = reads_pos.get(pos, 0) + 1
            reads_pos_score[pos] = reads_pos_score.get(pos, 0) + score

        for pos, base in enumerate(record.seq):
            if base in {"G", "C"}:
                gc_total += 1
                gc_per_position[pos] = gc_per_position.get(pos, 0) + 1
    avg_q_score_pos = {int(pos+1): v2/v1 for pos, (v1, v2) in enumerate(zip(reads_pos.values(), reads_pos_score.values()))}
    avg_gc_score_pos = {int(pos+1): round((v2/v1)*100, 2) for pos, (v1, v2) in enumerate(zip(reads_pos.values(), gc_per_position.values()))}

    return smallest_read, largest_read, read_lens, avg_q_score_read, gc_total, avg_q_score_pos, avg_gc_score_pos


def N50(read_lens):
    """
    Calculate N50 statistic for raw reads.
    """
    half_total_length = float(sum(read_lens) / 2)
    cumulative_length = 0
    for read in np.sort(read_lens)[::-1]:
        cumulative_length += read
        if cumulative_length >= half_total_length:
            return int(read)


def write_json(json_info, filename="QC_result.json"):
    """
    Writes all the statistics to a JSON file of youre choise.
    If not specified it will write to QC_report.json.
    """
    with open(file=f"results/QC{filename}", mode='w', encoding='utf-8') as f:
        json.dump(json_info, f, ensure_ascii=False, indent=4)


def run_analysis(fastq_file):
    """
    Collects all the statistics from the process_reads() and N50()
    and parses them to write_json().
    """
    (smallest_read, largest_read, read_lens, avg_q_score_read,
    gc_total, avg_q_score_pos, avg_gc_score_pos) = process_reads(fastq_file)
    total_bases = sum(read_lens)
    json_info = {
        "Smallest read": smallest_read,
        "Largest read": largest_read,
        "Mean read length": round(total_bases / len(read_lens), 2),
        "Total reads": len(read_lens),
        "Average GC total": round((gc_total / total_bases) * 100, 2),
        "Q-score every read": avg_q_score_read,
        "Average Q-score per position": avg_q_score_pos,
        "Average GC% per position": avg_gc_score_pos,
        "N50": N50(read_lens)
    }
    write_json(json_info)


if __name__ == "__main__":
    with gzip.open("data/reads/SRR8137614.fastq.gz", "rt") as fastq_file:
        run_analysis(fastq_file)
