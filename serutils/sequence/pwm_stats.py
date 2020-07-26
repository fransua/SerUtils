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

    h, mu  = zip(*sorted([(1 - abs(C1[i] - C2[i]), max(C1[i], C2[i])) 
                          for i in range(4)], reverse=True))
    mu = list(mu)

    x = newton(_func2, x0=-1, args=(mu,))

    muf = mu[0]
    for i in range(1, 4):
        muf = muf + mu[i] + x * muf * mu[i]
        mu[i] = muf
    
    sim = max(min(mu[i], h[i]) for i in range(4))
    return sim


def _func2(x, mu):
    return reduce(mul, (1 + v * x for v in mu)) - (1 + x)


def get_fisim(motif1, motif2, count_unaligned=False):
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
    :param False count_unaligned: if True, it will fill gaps in the shortest
       motif with equiprobable base frequencies. This penalizes motifs of
       different sizes.
    
    :return: FiSim score
    """
    motif1 = list(zip(*motif1))
    motif2 = list(zip(*motif2))

    if len(motif1) > len(motif2):
        motif1, motif2 = motif2, motif1
    diff = len(motif2) - len(motif1)

    all_fisim = []
    for i in range(diff + 1):
        this_fisim = 0
        for j in range(len(motif1)):
            C1 = motif1[j]
            C2 = motif2[j+i]
            this_fisim += _single_fisim(C1, C2)
        all_fisim.append((this_fisim, i))

    if not count_unaligned:
        return max(all_fisim)[0] / len(motif1)

    i = max(all_fisim)[1]
    this_fisim = 0
    for j in range(len(motif2)):
        if j < i or j >= i + len(motif1):
            C1 = [0.25] * 4
        else:
            C1 = motif1[j - i]
        C2 = motif2[j]
        this_fisim += _single_fisim(C1, C2)

    return this_fisim / len(motif2)
