from matplotlib import pyplot as plt
from matplotlib.patches import Ellipse
import numpy as np
from scipy.optimize import dual_annealing

def intersection_area(r1, r2, d):
    if d > r1 + r2:
        return 0
    r12 = r1**2
    r22 = r2**2

    if d <= abs(r2 - r1):
        return np.pi * min(r12, r22)
    
    d1 = (r12 - r22 + d**2) / (d * 2)
    d2 = d - d1
    
    A = (r12 * np.arccos(d1 / r1) - d1 * (r12 - d1**2)**0.5 + 
         r22 * np.arccos(d2 / r2) - d2 * (r22 - d2**2)**0.5)
    
    return A

def optimal_intersection(d, r1, r2, i):
    tmp_i = intersection_area(r1, r2, d)
    return abs(i - tmp_i)

def get_distance(r1, r2, i):
    return dual_annealing(optimal_intersection, args=(r1, r2, i,),
                           bounds=[(abs(r1 - r2), (r1 + r2))]).x


def scatter_venns(intersections, unions, sizes1, sizes2, xcoords,
                  ylim=(0, 1), factor=1, colors=None, axe=None, 
                  metric='Jaccard-index', title=None, scale_range=None):
    """
    plot overlap between datasets of size N, but each comparison
    is represented by a Venn-diagram.
    
    :param intersections: list of N comparisons with the number of overlaping elements.
    :param unions: list of N comparisons with total number of elements in the union.
    :param sizes1: list of N sizes of the first element in the comparisons.
    :param sizes2: list of N sizes of the second element in the comparisons.
    :param xcoords: list of X coordinates on to place the comparisons.
    :param (0, 1) ylim: Y limits
    :param 1 factor: scale factor for the size of the circles. Area of circles 
       are equal to the size of the element divided by twice the maximum size observed.
    :param None colors: list of pair of values to color each element in each comparison
    :param None axe: matplolib Axe object
    :param 'Jaccard-index' metric: metric to be used for the Y coordinate 
       (Jaccard-index or Overlap-coefficient)
    :param None title: plot title
    :param None scale_range: list of values to show examples of circle sizes in the legend.
       By default percentiles 1, 50 and 99 of the size distribution will be used.
    """
    if axe is None:
        fig = plt.figure(figsize=(8, 8), facecolor='w')
        axe1 = fig.add_subplot()
        axe2 = fig.add_subplot()
    else:
        axe1 = axe
        axe2 = axe.get_figure().add_subplot()

    if metric == 'Jaccard-index':
        metricf = lambda i, u, s1, s2: i / u
    elif metric == 'Overlap-coefficient':
        metricf = lambda i, u, s1, s2: i / min(s1, s2)
    else:
        raise NotImplementedError('WRONG metric: use "Jaccard-index" or "Overlap-coefficient"')

    if colors is None:
        colors = [[None, None]] * len(xcoords)

    factor = factor * max(np.max(sizes1), np.max(sizes2)) * 100
    
    xmin = min(xcoords)
    xmax = max(xcoords)
    xdiff = xmax - xmin
    xmin = xmin - xdiff * 0.15
    xmax = xmax + xdiff * 0.15
    xdiff = xmax - xmin
    ydiff = ylim[1] - ylim[0]
    
    Xratio = xdiff / ydiff
    
    for s1s2, (i, u, s1, s2, x, color) in enumerate(zip(
        intersections, unions, sizes1, sizes2, xcoords, colors)):
        ## stats/convertions
        # size is circle area, radius is thus square root of size over pi 
        r1 = (s1 / factor / np.pi)**0.5
        r2 = (s2 / factor / np.pi)**0.5
        # compute distance for venn-diagram
        d = get_distance(r1, r2, i / factor) / 2 * Xratio
        # get metric for Y coordinate (Jaccard-index by default)
        ii = metricf(i, u, s1, s2)
        
        # draw
        c = Ellipse((x - d, ii), r1 * Xratio * 2, r1 * 2, facecolor='none', 
                    edgecolor='k', linewidth=1, alpha=1, zorder=100)
        axe1.add_artist(c)
        c = Ellipse((x - d, ii), r1 * Xratio * 2, r1 * 2, 
                    facecolor=colors[s1s2][0], linewidth=1, alpha=0.5, zorder=100)
        axe1.add_artist(c)
        c = Ellipse((x + d, ii), r2 * Xratio * 2, r2 * 2, facecolor='none', 
                    edgecolor='k', linewidth=1, alpha=1, zorder=100)
        axe1.add_artist(c)
        c = Ellipse((x + d, ii), r2 * Xratio * 2, r2 * 2, 
                    facecolor=colors[s1s2][1], linewidth=1, alpha=0.5, zorder=100)
        axe1.add_artist(c)

    axe1_size = 0.6
    axe1.set_position((0.1, 0.1, axe1_size, axe1_size))
    axe1.set_title(title)
    # size legend
    size = xmax - xmin
    axe2_height = 0.17
    axe_ratio = axe1_size / axe2_height
    axe2.set_position((0.71, 0.05, 0.17, axe2_height))
    for n, s in enumerate(scale_range if scale_range else np.percentile(
        list(sizes1) + list(sizes2), (1, 50, 99))):
        y = np.linspace(ylim[0], ylim[1], 15)[n * 5 + 1]
        r = (axe1_size / axe2_height) * (s / factor / np.pi)**0.5
        e = Ellipse((xmin + size / 5, y), r * Xratio * 2, r * 2,
                    facecolor='none', edgecolor='k', linewidth=1, clip_on=False)
        axe2.text(0.9, y, f'{s // 1000}K', va='center')
        axe2.add_artist(e)
    axe2.set_xlim(-size / 2, size / 2)
    axe2.set_ylim(ylim)
    axe2.set_axis_off()

    axe2.set_title('    Size scale', ha='center')

    axe1.set_xlim(xmin, xmax)
    axe1.set_ylim(ylim)
    axe1.set_ylabel(metric)
    axe1.grid()
