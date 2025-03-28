from collections import Counter
import difflib

class HaplotypeAligner:
    def __init__(self, haplotype1, haplotype2, hap1_is_negative=False, hap2_is_negative=False):
        self.haplotype1 = haplotype1[::-1] if hap1_is_negative else haplotype1
        self.haplotype2 = haplotype2[::-1] if hap2_is_negative else haplotype2

        self.hap1_segment_count = Counter([segment[0] for segment in self.haplotype1])
        self.hap2_segment_count = Counter([segment[0] for segment in self.haplotype2])

        self.hap1_seen_segments = {}
        self.hap2_seen_segments = {}

        self.hap1_unaligned = []
        self.hap2_unaligned = []

        self.alignment_result = []

        self.matcher = difflib.SequenceMatcher(None, self.haplotype1, self.haplotype2)


    def get_duplicate_segment_name(self, segment_name, seen_counter):
        """Adds numbering to duplicate segment names."""
        seen_counter[segment_name] = seen_counter.get(segment_name, 0) + 1
        return f"{segment_name} ({seen_counter[segment_name]})"


    def flush_unaligned_segments(self):
        """Places unaligned segments in the correct order."""
        while self.hap1_unaligned or self.hap2_unaligned:
            if not self.hap1_unaligned:
                self.alignment_result.append((None, self.hap2_unaligned.pop(0)[0]))
            elif not self.hap2_unaligned:
                self.alignment_result.append((self.hap1_unaligned.pop(0)[0], None))
            else:
                hap1_segment, hap1_start, hap1_end = self.hap1_unaligned[0]
                hap2_segment, hap2_start, hap2_end = self.hap2_unaligned[0]

                if hap1_end < hap2_start:
                    self.alignment_result.append((self.hap1_unaligned.pop(0)[0], None))
                elif hap2_end < hap1_start:
                    self.alignment_result.append((None, self.hap2_unaligned.pop(0)[0]))
                else:
                    self.alignment_result.append((self.hap1_unaligned.pop(0)[0], self.hap2_unaligned.pop(0)[0]))


    def align(self):
        """Aligns segments from two haplotypes."""
        for operation_type, hap1_start, hap1_end, hap2_start, hap2_end in self.matcher.get_opcodes():
            if operation_type == "equal":
                for idx1, idx2 in zip(range(hap1_start, hap1_end), range(hap2_start, hap2_end)):
                    self.process_segment_pair(idx1, idx2)
            elif operation_type == "delete":
                self.process_deletion(hap1_start, hap1_end)
            elif operation_type == "insert":
                self.process_insertion(hap2_start, hap2_end)
            elif operation_type == "replace":
                self.process_replacement(hap1_start, hap1_end, hap2_start, hap2_end)

        self.flush_unaligned_segments()
        return self.alignment_result


    def process_segment_pair(self, hap1_index, hap2_index):
        """Processes a pair of aligned segments."""
        hap1_segment, hap1_start, hap1_end = self.haplotype1[hap1_index]
        hap2_segment, hap2_start, hap2_end = self.haplotype2[hap2_index]

        # Handle duplicate segment names
        if self.hap1_segment_count[hap1_segment] > 1:
            hap1_segment = self.get_duplicate_segment_name(hap1_segment, self.hap1_seen_segments)
        if self.hap2_segment_count[hap2_segment] > 1:
            hap2_segment = self.get_duplicate_segment_name(hap2_segment, self.hap2_seen_segments)

        if hap1_end < hap2_start or hap2_end < hap1_start:
            self.hap1_unaligned.append((hap1_segment, hap1_start, hap1_end))
            self.hap2_unaligned.append((hap2_segment, hap2_start, hap2_end))
        else:
            self.flush_unaligned_segments()
            self.alignment_result.append((hap1_segment, hap2_segment))


    def process_deletion(self, hap1_start, hap1_end):
        """Processes segments present in haplotype1 but missing in haplotype2."""
        for hap1_index in range(hap1_start, hap1_end):
            hap1_segment, hap1_start, hap1_end = self.haplotype1[hap1_index]
            if self.hap1_segment_count[hap1_segment] > 1:
                hap1_segment = self.get_duplicate_segment_name(hap1_segment, self.hap1_seen_segments)
            self.hap1_unaligned.append((hap1_segment, hap1_start, hap1_end))


    def process_insertion(self, hap2_start, hap2_end):
        """Processes segments present in haplotype2 but missing in haplotype1."""
        for hap2_index in range(hap2_start, hap2_end):
            hap2_segment, hap2_start, hap2_end = self.haplotype2[hap2_index]
            if self.hap2_segment_count[hap2_segment] > 1:
                hap2_segment = self.get_duplicate_segment_name(hap2_segment, self.hap2_seen_segments)
            self.hap2_unaligned.append((hap2_segment, hap2_start, hap2_end))


    def process_replacement(self, hap1_start, hap1_end, hap2_start, hap2_end):
        """Processes segments that are replaced between haplotype1 and haplotype2."""
        min_length = min(hap1_end - hap1_start, hap2_end - hap2_start)
        for offset in range(min_length):
            self.process_segment_pair(hap1_start + offset, hap2_start + offset)

        if hap1_end > hap1_start + min_length:
            self.process_deletion(hap1_start + min_length, hap1_end)
        elif hap2_end > hap2_start + min_length:
            self.process_insertion(hap2_start + min_length, hap2_end)
