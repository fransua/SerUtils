
from collections import OrderedDict


def parse_fasta(f_name, verbose=True):
    """
    Parse one fasta.
    :param f_name: single path

    :returns: a sorted dictionary with chromosome names as keys, and sequences
       as values (sequence in upper case)
    """
    genome_seq = OrderedDict()
    for line in open(f_name):
        if line.startswith('>'):
            header = line[1:].split()[0]
            genome_seq[header] = ''
            if verbose:
                print('Parsing %s' % (header))
            continue
        genome_seq[header] += line.rstrip()
    return genome_seq
