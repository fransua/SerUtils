"""
"""
import numpy as np
from scipy.optimize import newton
from operator import mul
from functools import reduce


def _single_fisim(C1, C2):
    """
    """
    if C1 == C2:
        return 1

    size = len(C1)

    h, mu  = zip(*sorted([(1 - abs(C1[i] - C2[i]), max(C1[i], C2[i]))
                          for i in range(size)], reverse=True))
    mu = list(mu)

    x = newton(_func2, x0=-1, args=(mu,))

    muf = mu[0]
    for i in range(1, size):
        muf = muf + mu[i] + x * muf * mu[i]
        mu[i] = muf
    
    sim = max(min(mu[i], h[i]) for i in range(size))
    return sim


def _func2(x, mu):
    return reduce(mul, (1 + v * x for v in mu)) - (1 + x)


def get_fisim(motif1, motif2, check_reverse=True, count_unaligned=False):
    """
    Get the Fuzzy Integral SIMilarity measure (Garcia et al. 2009) for a pair
    of DNA binding motifs. The function computes the FiSim at all possible
    positions and peaks the best.

        Garcia, Fernando, et al. 
        "FISim: a new similarity measure between transcription factor binding
        sites based on the fuzzy integral." 
        BMC bioinformatics 10.1 (2009): 224.

    :param motif1: PWM of the first motif (matrix of 4 lines and as many
       column as motif length)
    :param motif2: PWM of the first motif (matrix of 4 lines and as many
       column as motif length)
    :param True check_reverse: if True, it will try to align shortest motif
       in reverse orientation also.
    :param False count_unaligned: if True, it will fill gaps in the shortest
       motif with equiprobable base frequencies. This penalizes motifs of
       different sizes.
    
    :return: FiSim score between 0 and 1 (1:identical; 0 different), and index
       (negative index if best matching was found reversing shortest motif)
    """
    motif1 = list(zip(*motif1))
    motif2 = list(zip(*motif2))

    if len(motif1) > len(motif2):
        motif1, motif2 = motif2, motif1
    diff = len(motif2) - len(motif1)

    all_fisim = []
    computed = {}  # for memoization
    strands = [1, -1] if check_reverse else [1]
    # once, both forward, once reversing the shortest motif
    max_fisim = 0
    for strand in strands:
        # slide short motif along long motig
        for i in range(diff + 1):
            this_fisim = 0
            # iterate over columns in the alignment
            for j in range(len(motif1)):
                C1 = tuple(motif1[::strand][j])
                C2 = tuple(motif2[j+i])
                try:
                    this_fisim += computed[C1, C2]
                except KeyError:
                    fisim = _single_fisim(C1, C2)
                    computed[C1, C2] = fisim
                    this_fisim += fisim
                # stop if we already lost the game
                if this_fisim + (len(motif1) - j - 1) < max_fisim:
                    break
            all_fisim.append((this_fisim, strand * i))
            max_fisim = max(all_fisim)[0]

    if not count_unaligned:
        fisim, i = max(all_fisim)
        return  fisim / len(motif1), i

    i = abs(max(all_fisim)[1])
    this_fisim = 0
    for j in range(len(motif2)):
        if j < i or j >= i + len(motif1):
            C1 = tuple([0.25] * len(motif1[0]))
        else:
            C1 = tuple(motif1[j - i])
        C2 = tuple(motif2[j])
        try:
            this_fisim += computed[C1, C2]
        except KeyError:
            fisim = _single_fisim(C1, C2)
            computed[C1, C2] = fisim
            this_fisim += fisim

    return this_fisim / len(motif2), max(all_fisim)[1]
