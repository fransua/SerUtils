"""
"""

try:
    from itertools import zip_longest
except ImportError:
    from itertools import izip_longest as zip_longest
from random import random, seed
from bisect import bisect_right as bisect


def generate_random_genomic_positions(coordinates=None, sequences=None,
                                      num=1, seed_num=1):
    seed(seed_num)
    if not coordinates:
        if sequences:
            coordinates = dict((c, len(sequences[c])) for c in sequences)
        else:
            raise Exception('ERROR: at least one of coordinates or sequences should be given!!')

    total_length = sum(v for v in coordinates.values())

    positions = [0]
    chroms = []
    total = 0
    for c in sorted(coordinates.keys(), key=lambda x: coordinates[x]):
        total += coordinates[c]
        positions.append(total)
        chroms.append(c)

    if sequences:
        count = 0
        while count < num:
            pos = int(random() * total_length)
            idx = bisect(positions, pos) - 1
            c, pos = chroms[idx], pos - positions[idx]
            if sequences[c][pos] == 'N':
                continue
            yield (c, pos + 1)
            count += 1
    else:
        for _ in range(num):
            pos = int(random() * total_length)
            idx = bisect(positions, pos) - 1
            yield (chroms[idx], pos - positions[idx] + 1)


def nucleotide_content(seqs, letters=('A', 'T', 'G', 'C', 'N'), remove_gap_only=True, min_count=1):
    """
    :params None seqs: list of sequences
    :param 'A', 'T', 'G', 'C', 'N' letters: nucleotides to consider

    :returns: list of dictionaries of nucleotide proportions per sequence.
    """
    seq_count = []
    col_count = []
    lencol = float(len(seqs))
    for col in zip_longest(*seqs, fillvalue='N'):
        nchars = lencol - col.count('-')
        if nchars < min_count:
            if not remove_gap_only:
                seq_count.append(dict([(let, 0) for let in letters]))
                col_count.append(nchars)
            continue
        col_count.append(nchars)
        seq_count.append(dict([(let, col.count(let) / nchars) for let in letters]))
    return seq_count, col_count
