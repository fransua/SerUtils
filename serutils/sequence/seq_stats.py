from itertools import izip_longest

import numpy as np
from matplotlib import pyplot as plt


def seq_logo(seqs, title='', axe=None, quality=1, plot=True, ylim=(0, 2), lwmodif=1.25,
             errorbar=True, savefig=None):
    """
    Draw a sequence motif according to an input list of sequences,
    https://en.wikipedia.org/wiki/Sequence_logo

    :params seqs: list of sequences
    :params '' title: title of the plot
    :params None axe: matplotlib axe object
    :params 1 quality: quality modifier
    :params None savefig: if set save the plot to the file at the specified path
       (the format of the file "pdf", "png"... will be infered from the file extension)
    :params 1.25 lwmodif: modifier for linewith of the drown letters.
    :params True plot: draw plot, otherwise only returns proportions.
    :params True errorbar: draw error bar, that are equal to the correction for
       small samples

    :returns: a table (list of list) with the proportions of each nucleotide at
       each position
    """
    xoffset = 0.003 * quality
    yoffset = 0.030 * quality
    lwmodif = 1.25
    letterw = 0.75

    qlmodif = lwmodif * quality / 5.
    halfw   = letterw / 2.

    letters = ['A', 'T', 'G', 'C', 'N']

    nseq = float(len(seqs))
    seq_count = [dict([(let, col.count(let) / nseq) for let in letters])
                 for col in izip_longest(*seqs, fillvalue='N')]

    if not plot:
        return seq_count

    def lw(val):
        return ((val + np.log(val) + 7) * qlmodif) if val else 0.

    def A(x, y, h):
        if yoffset < h:
            y += yoffset
            h -= yoffset
        axe.plot([x, x + halfw, x + letterw ], [y, y + h, y], color='green', lw=lw(h))
        axe.plot([x + letterw * 1 / 4, x + letterw * 3 / 4], [y + h  * 1 / 3] * 2,
                 color='green', lw=lw(h))

    def T(x, y, h):
        if yoffset < h:
            y += yoffset
            h -= yoffset
        axe.plot([x + halfw ] * 2, [y, y + h] , color='red', lw=lw(h))
        axe.plot([x, x + letterw]  , [y + h] * 2, color='red', lw=lw(h))

    Gcircle = np.linspace(np.pi * 0.25, np.pi * 1.9, 20 * quality)
    xG = np.cos(Gcircle)
    yG = np.sin(Gcircle)
    def G(x, y, h):
        if yoffset < h:
            y += yoffset
            h -= yoffset
        halfh = h / 2.
        xGw = x + xG * halfw + halfw
        yGh = y + yG * halfh + halfh
        axe.plot(xGw, yGh, color='orange', lw=lw(h))
        axe.plot([xGw[-1]          , xGw[-1]           , xGw[-1] - letterw * 0.4 ],
                 [yGh[-1] - h * 0.3, yGh[-1] + h * 0.05, yGh[-1] + h       * 0.05],
                color='orange', lw=lw(h))

    Ccircle = np.linspace(np.pi * 0.2, np.pi * 1.8, 20 * quality)
    xC = np.cos(Ccircle)
    yC = np.sin(Ccircle)
    def C(x, y, h):
        if yoffset < h:
            y += yoffset
            h -= yoffset
        axe.plot(x + xC * halfw + halfw,
                 y + yC * h / 2. + h / 2., color='blue', lw=lw(h))

    def N(x, y, h):
        if yoffset < h:
            y += yoffset
            h -= yoffset
        axe.plot([x, x, x + letterw, x + letterw], [y, y + h, y, y + h], color='grey', lw=lw(h))

    if not axe:
        fig = plt.figure(figsize=(len(seq_count) * 0.15 * quality, 1.75 * quality))
        axe = fig.add_subplot(111)

    draw_letter_functions = {'A':A, 'T': T, 'G': G, 'C': C, 'N': N}

    en = 1. / np.log(2) * (4 - 1.) / (2. * nseq)
    for npos, position in enumerate(seq_count):
        offset = 0
        ri = sum(fi * np.log2(fi) for fi in position.values() if fi) - en
        for letter, f in sorted(position.iteritems(), key=lambda x: x[1]):
            h = f * (2. + ri)
            x = npos + xoffset
            y = offset
            draw_letter_functions[letter](x, y, h)
            if errorbar:
                plt.plot([x + halfw, x + halfw], [y + h + en, max(y + h - en, 0)],
                         color='grey')
                plt.plot([x + halfw - 0.1, x + halfw + 0.1], [y + h + en, y + h + en],
                         color='grey')
                if y + h - en > 0:
                    plt.plot([x + halfw - 0.1, x + halfw + 0.1], [y + h - en] * 2,
                             color='grey')
            offset += h

    # axes
    axe.axis('off')
    axe.plot([-0.4, -0.1, -0.1, -0.4], [ylim[0], ylim[0], ylim[1], ylim[1]], color='k')
    for i in range(3):
        axe.text(-0.5, ylim[0] + (ylim[1] - ylim[0]) / 2. * i,
                 str(ylim[0] + (ylim[1] - ylim[0]) / 2. * i),
                 va='center', ha='right', size=4.5 * quality)
    axe.text(ylim[0] - 1.5, 1, 'bits', va='center', ha='right', rotation=90, size=7 * quality)
    axe.plot([-0.4, -0.1], [(ylim[1] - ylim[0]) / 2.] * 2, color='k')
    for i in xrange(1, npos + 2):
        axe.text(i - 0.5, ylim[0] - 0.05, str(i), va='top', ha='center',
                 size=4.5 * quality, rotation=90)
    _ = axe.set_xlim(-1, npos + 1.)
    _ = axe.set_ylim(ylim[0] - .5, ylim[1] + 3 * yoffset + 0.5)
    # title
    axe.set_title(title)
    if savefig:
        plt.savefig(savefig, format=savefig.split('.')[-1])
    return seq_count
