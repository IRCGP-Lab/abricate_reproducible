from collections import defaultdict
from dataclasses import dataclass
from typing import List
import argparse


@dataclass
class AbricateEntry:
    def __init__(self, parse_entry_line: str, delimiter = ","):
        self.line = parse_entry_line
        self._entryline_parse(parse_entry_line, delimiter)

    def _entryline_parse(self, parse_entry_line: str, delimiter: str):
        """
        0: #FILE, 1: SEQUENCE, 2: START, 3: END, 4: STRAND, 5: GENE, 6: COVERAGE, 7: COVERAGE_MAP,
        8: GAPS, 9 :%COVERAGE, 10 : %IDENTITY, 11: DATABASE, 12: ACCESSION, 13: PRODUCT, 14: RESISTANCE
        :param parse_entry_line:
        :param delimiter:
        :return:
        """
        entry = parse_entry_line.split(delimiter)
        self.contig = entry[1]
        self.start = int(entry[2])
        self.end = int(entry[3])
        self.coverage = float(entry[9]) * 100
        self.identity = float(entry[10]) * 100
        self.score = self.coverage * self.identity

    def __repr__(self):
        return self.line

    # Rich comparison methods
    # def __lt__(self, other):
    #     return self.start < other.start
    #
    # def __eq__(self, other):
    #     if not isinstance(other, self.__class__):
    #         return False
    #     return (self.start == other.start) and (self.end == other.end) and (self.contig == other.contig)

    @classmethod
    def _group_reads(cls, entries: List['abricate_entry']):
        sorted_entries = sorted(entries, key=lambda x: x.start)
        start = -1
        end = -1
        interval_bubbles = []
        bubble_group = []
        for entry in sorted_entries:
            if start <= entry.start < end:
                if entry.end > end:
                    end = entry.end
            elif start < entry.end <= end:
                if entry.start < start:
                    start = entry.start
            elif (start >= entry.start) and (end <= entry.end):
                start = entry.start
                end = entry.end
            else:
                if bubble_group:
                    interval_bubbles.append(bubble_group)
                    bubble_group = []
                start = entry.start
                end = entry.end
            bubble_group.append(entry)
        if bubble_group:
            interval_bubbles.append(bubble_group)
        results = []
        for bubble in interval_bubbles:
            sorted_bubble = sorted(bubble, key=lambda x: x.score, reverse=True)
            best_score = sorted_bubble[0].score
            for entry in sorted_bubble:
                if entry.score == best_score:
                    results.append(entry)
                elif entry.score > best_score:
                    raise ArithmeticError
        return results

    @classmethod
    def summarize(cls, results: str, output: str, delimiter=","):
        contig_wise = defaultdict(list)
        header = ""
        with open(results) as abr:
            for line in abr:
                if not line.startswith("#"):
                    contig_id = line.split(delimiter)[1]
                    strand = line.split(delimiter)[4]
                    id_strand = contig_id + strand
                    contig_wise[id_strand].append(AbricateEntry(line, delimiter=delimiter))
                else:
                    header = line
        result_report = []
        for contig in contig_wise:
            result_report.extend(cls._group_reads(contig_wise[contig]))
        with open(output, "w") as o:
            o.write(header)
            for i in result_report:
                o.write(i.line)


def main():
    parser = argparse.ArgumentParser(description="abricate parser")

    parser.add_argument('input', type=str, help='raw abricate result')
    parser.add_argument('output', type=str, help='processed file')
    parser.add_argument('delimiter', type=str, default=",", help='delimiter')

    args = parser.parse_args()

    AbricateEntry.summarize(results=args.input, output=args.output, delimiter=args.delimiter)


if __name__ == '__main__':
    main()

