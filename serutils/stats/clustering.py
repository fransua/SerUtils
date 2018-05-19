"""
"""

from itertools               import combinations

from scipy.cluster.hierarchy import fcluster


def calinski_harabasz(scores, clusters):
    """
    Implementation of the CH score [CalinskiHarabasz1974]_, that has shown to be
    one the most accurate way to compare clustering methods
    [MilliganCooper1985]_ [Tibshirani2001]_.

    The CH score is:

    .. math::

        CH(k) = \\frac{B(k) / (k-1)}{W(k)/(n - k)}

    Where :math:`B(k)` and :math:`W(k)` are the between and within cluster sums
    of squares, with :math:`k` clusters, and :math:`n` the total number of
    elements (models in this case).

    :param scores: a dict with, as keys, a tuple with a pair of models; and, as
       value, the distance between these models.
    :param clusters: a dict with, as key, the cluster number, and as value a
       list of models
    :param nmodels: total number of models

    :returns: the CH score
    """
    cluster_list = [c for c in clusters if len(clusters[c]) > 1]
    if len(cluster_list) <= 1:
        return 0.0
    nmodels = sum([len(clusters[c]) for c in cluster_list])

    between_cluster = (sum([sum([sum([scores[(md1, md2)]**2
                                      for md1 in clusters[cl1]])
                                 for md2 in clusters[cl2]])
                            / (len(clusters[cl1]) * len(clusters[cl2]))
                            for cl1, cl2 in combinations(cluster_list, 2)])
                       / ((len(cluster_list) - 1.0) / 2))

    within_cluster = (sum([sum([scores[(md1, md2)]**2
                                for md1, md2 in combinations(clusters[cls], 2)])
                           / (len(clusters[cls]) * (len(clusters[cls]) - 1.0) / 2)
                           for cls in cluster_list]))

    return ((between_cluster / (len(cluster_list) - 1))
            /
            (within_cluster / (nmodels - len(cluster_list))))


def find_best_clusters(hcluster, score_matrix, leaves):
    # score each possible cut in hierarchical clustering
    solutions = {}
    for k in hcluster[:, 2]:
        clusters = {}
        [clusters.setdefault(j, []).append(i) for i, j in
         enumerate(fcluster(hcluster, k, criterion='distance'))]
        solutions[k] = {'out': clusters}
        solutions[k]['score'] = calinski_harabasz(score_matrix, clusters)
    # take best cluster according to calinski_harabasz score
    clusters = [solutions[s] for s in sorted(
        solutions, key=lambda x: solutions[x]['score'])
        if solutions[s]['score'] > 0][-1]['out']
    # sort clusters, the more populated, the first.
    clusters = dict([(i + 1, j) for i, j in
                     enumerate(sorted(clusters.values(),
                                      key=len, reverse=True))])
    new_clusters = {}
    leaf_cluster = {}
    for cluster in clusters:
        new_clusters[cluster] = []
        for leaf in clusters[cluster]:
            leaf_cluster[leaf] = cluster
            new_clusters[cluster].append(str(leaf))
        # new_clusters[cluster].sort(
        #     key=lambda x: self[str(x)]['objfun'])

    dads = {}
    i = max(clust_count)
    for a, b, _, _ in z:
        i += 1
        clust_count[i] = clust_count[a + 1] + clust_count[b + 1]
        dads[a + 1] = i
        dads[b + 1] = i

    d = augmented_dendrogram(clust_count, dads, objfun, color,
                             axe, savefig, z, **kwargs)
    return d


def main():
    from scipy.cluster.hierarchy import linkage

    
