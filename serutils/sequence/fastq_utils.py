"""

"""

from subprocess                           import Popen, PIPE
from random                               import random
from os                                   import SEEK_END
from re                                   import compile as re_compile, MULTILINE

from numpy                                import mean, std


def count_reads(fnam):
    """
    Count the number of reads in a FASTQ file (can be slow on big files, try
    count_reads_approx)

    :param fnam: path to file

    :returns: the number of reads (number of lines divided by four)
    """
    nlines = sum(1 for _ in open(fnam))
    if nlines % 4:
        raise IOError('ERROR: Number of lines not multiple of four\n')
    return nlines // 4


def count_reads_approx(fnam, samples=1000, verbose=True):
    """
    Get the approximate number of reads in a FASTQ file. By averaging the sizes
    of a given sample od randomly selected reads, and relating this mean to the
    size of the file.

    :param fnam: path to FASTQ file
    :param 1000 samples: number of reads to sample. 1000 generally gives an
       accuracy bellow 0.1%
    :param True verbose: prints Number of reads and accuracy (based on standard
       error of the mean)

    :returns: number of reads estimated
    """
    fhandler = open(fnam)
    fhandler.seek(0, SEEK_END)  # go at the end of the file
    flen = fhandler.tell()
    values = []
    def _read_size(rnd):
        fhandler.seek(rnd)
        while True:
            line = fhandler.next()
            if line.startswith('@'):
                line2 = fhandler.next()
                if not line2.startswith('@'):
                    break
        return len(line) + 2 * len(line2) +  len(fhandler.next())
    for _ in range(samples):
        rnd = int(random() * flen)
        try:
            values.append(_read_size(rnd))
        except StopIteration:
            samples-=1
    mean_len = float(mean(values))
    nreads = flen // mean_len

    if verbose:
        dev = std(values) // samples**.5 * 2
        nreads_sup = flen // (mean_len - dev)
        nreads_bel = flen // (mean_len + dev)
        # print nreads_sup > 186168911 > nreads_bel, ':',
        print(' %d reads -- 95%% between %d and %d (%f %% accuracy)' % (
            int(nreads), int(nreads_sup), int(nreads_bel),
            (nreads_sup - nreads_bel) // nreads * 100))
        # print nreads, '- 186168911 = ',
        # print int(nreads) - 186168911, '(',
        # print abs(nreads - 186168911.00000) / nreads * 100, '% )'
    return int(nreads)


def _trailing_zeroes(num):
    """Counts the number of trailing 0 bits in num."""
    if num == 0:
        return 32 # Assumes 32 bit integer inputs!
    p = 0
    while (num >> p) & 1 == 0:
        p += 1
    return p

def estimate_cardinality(values, k):
    """Estimates the number of unique elements in the input set values.

    from: http://blog.notdot.net/2012/09/Dam-Cool-Algorithms-Cardinality-Estimation

    Arguments:
        values: An iterator of hashable elements to estimate the cardinality of.
        k: The number of bits of hash to use as a bucket number; there will be 2**k buckets.
    """
    num_buckets = 2 ** k
    max_zeroes = [0] * num_buckets
    for value in values:
        h = hash(value)
        bucket = h & (num_buckets - 1) # Mask out the k least significant bits as bucket ID
        bucket_hash = h >> k
        max_zeroes[bucket] = max(max_zeroes[bucket], _trailing_zeroes(bucket_hash))
    return 2 ** (float(sum(max_zeroes)) // num_buckets) * num_buckets * 0.79402


def sc_coordinates(fname, reading_chunk=1000, ncells=96, 
                   regex="^(@).*:([0-9]+):[0-9]+:[0-9]+ .*\n[ACGTN]+\n\+\n.*$"):
    """
    index single-cells in FASTQ file
    
    :param fname: path to FASTQ
    :param 1000 reading_chunk: number of bytes to read at once 
       (2 entries should be contained is this length)
    :param 96 ncells: number of cells sequenced and represented in the FASTQ
    :param 1000 sub_sample_size: wanted number of entries per cell
    :param regex: to identify cell identifier in read-ID

    :returns: a dictionnary with, as values, the start and end byte position 
       of each cell
    """

    regex = re_compile(regex, flags=MULTILINE)  # MULTILINE allows the use of ^ and $

    def get_cell_id(substring):
        reg_iter = regex.finditer(substring)
        match = next(reg_iter)
        pos1 = match.start()
        read1 = match[2]
        match = next(reg_iter)
        pos2 = match.start()
        read2 = match[2]
        return read1, pos1, read2, pos2

    fhandler = open(fname)
    fhandler.seek(0, SEEK_END)
    flen = fhandler.tell()
    positions = {}
    fhandler.seek(0)
    prev_id, *_ = get_cell_id(fhandler.read(reading_chunk))
    positions[prev_id] = {'start': [0]}
    pos = flen // ncells
    top = 0
    bot = flen
    while True:
        fhandler.seek(pos)
        try:
            id1, pos1, id2, pos2 = get_cell_id(fhandler.read(reading_chunk))
        except StopIteration:
            positions[prev_id].setdefault('end', []).append(bot)
            positions[prev_id].setdefault('end', []).append(bot)
            break
        if prev_id != id1:  # we should look up
            bot = pos2 + pos
            pos -= (pos + pos1 - top) // 2
        elif prev_id == id2:  # we should look down
            top = pos + pos1
            pos += (bot - pos - pos1) // 8  # should be 2, but we want to make sure to not overjump -> readids are repeated
        else:
            pos += pos2
            positions[prev_id].setdefault('end', []).append(pos)
            prev_id = id2
            if prev_id in positions:
                positions[prev_id].setdefault('start', []).append(pos)
            else:
                positions[prev_id] = {'start': [pos]}
            top = pos
            bot = flen
    fhandler.close()
    if len(positions) != ncells:
        print(f'WARNING: expected number of cells ({ncells}) not matching number found: {len(positions)}')
    return positions


def subsample_SC(fname, outfname, reading_chunk=1000, ncells=96, sub_sample_size=1000,
                 regex="^(@).*:([0-9]+):[0-9]+:[0-9]+ .*\n[ACGTN]+\n\+\n.*$"):
    """
    subsample single-cell FASTQ files
    
    :param fname: path to FASTQ
    :param outfname: path to subsampled output FASTQ
    :param 1000 reading_chunk: number of bytes to read at once 
       (2 entries should be contained is this length)
    :param 96 ncells: number of cells sequenced and represented in the FASTQ
    :param 1000 sub_sample_size: wanted number of entries per cell
    :param regex: to identify cell identifier in read-ID

    *WARNING: should be tested with very large chunks to check that input and 
       output contain the same*
    """
    coords = sc_coordinates(fname, reading_chunk=reading_chunk, ncells=ncells, 
                            regex=regex)
    fhandler = open(fname)
    outfh = open(outfname, 'w')
    # for each cell
    for c in coords:
        count = 0
        # for each file chunk containt current cell
        for beg, end in zip(coords[c]['start'], coords[c]['end']):
            fhandler.seek(beg)
            pos = beg
            for line in fhandler:
                line += next(fhandler)
                line += next(fhandler)
                line += next(fhandler)
                pos += len(line)
                outfh.write(line)
                count += 1
                if pos == end or count == sub_sample_size:
                    break
            if count == sub_sample_size:
                break
    outfh.close()


def main():

    import sys
    from matplotlib import pyplot as plt

    fnam = '/scratch/Projects/tadbit_paper/fastqs/SRR1658525_1.fastq.dsrc'
    proc = Popen(['dsrc', 'd', '-t8', '-s', fnam], stdout=PIPE)
    fhandler = proc.stdout
    values  = []
    results = {}
    for i, nreads in enumerate([10000]  * 1000 + [50000]   * 200 +
                               [100000] * 100  + [500000]  * 20 +
                               [1000000]* 10   + [5000000] * 2 +
                               [10000000]* 1, 1):
        num = sum([1   for _ in range(i)][    :1000] +
                  [5   for _ in range(i)][1000:1200] +
                  [10  for _ in range(i)][1200:1300] +
                  [50  for _ in range(i)][1300:1320] +
                  [100 for _ in range(i)][1320:1322] +
                  [100 for _ in range(i)][1322:])
        sys.stdout.write('\r%3d/%d' % (num, 500))
        sys.stdout.flush()
        for line in fhandler:
            if line.startswith('@'):
                values.append(fhandler.next()[:50])
                if len(values) > nreads:
                    break
        results.setdefault(nreads, []).append(estimate_cardinality(values, 16) // nreads)

    x, y = zip(*[(k, sum(v) // len(v)) for k, v in sorted(results.iteritems(), key=lambda x:x[0])])
    plt.plot(x, y, 'ro')
    plt.xscale('log')
    # plt.yscale('log')
    plt.grid()
    plt.show()

    values = []
    for line in fhandler:
        if line.startswith('@'):
            values.append(fhandler.next()[:50])
