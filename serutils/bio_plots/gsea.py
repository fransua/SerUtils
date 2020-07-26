#!/usr/bin/env python

from cPickle import load
from collections import OrderedDict
from glob import glob
from IPython.display import display
import os

from pysam import AlignmentFile

import numpy as np

from matplotlib import pyplot as plt

from scipy.cluster.hierarchy import linkage, dendrogram, cut_tree
from scipy.spatial.distance import pdist, squareform

from goatools.godag_plot    import plot_gos, plot_results, plot_goid2goobj
from goatools.obo_parser    import GODag
from goatools.go_enrichment import GOEnrichmentStudy
from goatools.base          import download_go_basic_obo
from goatools.associations  import read_ncbi_gene2go
from goatools.associations  import read_gaf
from goatools.test_data.genes_NCBI_9606_ProteinCoding import GENEID2NT as GeneID2nt
from goatools.semantic import semantic_similarity, semantic_distance, lin_sim
from goatools.semantic import resnik_sim, lin_sim
from goatools.semantic import TermCounts
import goatools

goatools.__version__

cutlvl   = 12
strategy = "Lpv0.1"# "pv0.1" # 'switch'
#strategy = "switch"# "pv0.1" # 'switch'
# strategy = "ipv0.5"# "pv0.1" # 'switch'
goset    = 'goslim' # 'goslim'
goset    = 'allGOs' # 'goslim'
from_anc = 'nop' # primate, anc or nop


for spe in species.keys():
    short = spe.split()[0][0] + spe.split()[1][:3]
    # species[short] = species[spe]
    species[spe]['short'] = short


bam_path = '/mnt/cnag_scratch/Projects/GEVO_3D/results/{}/03_filtered_reads/intersection_*.bam'
for spe in species:
    bam = AlignmentFile(glob(bam_path.format(species[spe]['short']))[0])
    species[spe]['chrom_lens'] = dict(zip([c.upper() for c in bam.references], bam.lengths))


reso = 100000


out = open('/data/Projects/GEVO_3D/results/compartments/compare_AB/compareAB_100kb.pickle')
results = load(out)
out.close()


genes = {}
fh = open('/scratch/db/Genomes/gene_positions_Ensembl_v90/Hsap_ENSv90_GRCh38.p10.tsv')
fh.next()
for line in fh:
    try:
        _, gene_name, c, b, e, s = line.split('\t')
    except ValueError:
        print (line)
        continue
    genes.setdefault((c, (int(b) + int(e)) / 2 / reso), []).append(gene_name)
for k in genes:
    genes[k] = list(set(genes[k]))

out = open('/scratch/db/Genomes/gene_positions_Ensembl_v90/Hsap_100kb_ENSv90_GRCh38.p10.tsv', 'w')
for k in genes:
    out.write('%s\t%s\n' % ('%s:%d' % (k), ';'.join(genes[k])))
out.close()


lost_genes = load(open('EnrichR_lost_genes_100kb.pickle'))

for spe in species:
    if not spe == 'Homo sapiens':
        lost_genes[spe]  = set(lost_genes[spe])

get_ipython().system(' wget http://www.geneontology.org/ontology/subsets/goslim_generic.obo')


obo_fname = download_go_basic_obo()


from goatools.base import download_ncbi_associations


gene2go = download_ncbi_associations()


if goset=='goslim':
    obodag = GODag("goslim_generic.obo")
else:
    obodag = GODag("go-basic.obo")


geneid2gos = read_ncbi_gene2go("gene2go", taxids=[9606])


levels = [r.level for go, r in obodag.items()]
[(i, levels.count(i)) for i in range(1, 12)]


bad_go = []
for go, r in obodag.iteritems():
    if r.level > cutlvl:
        bad_go.append(go)
bad_go = set(bad_go)


len(bad_go)


for go, r in obodag.items():
    nps = set()
    for p in r._parents:
        if p in bad_go:
            nps |= set([p2 for p2 in obodag[p]._parents if not p2 in bad_go])
        else:
            nps.add(p)
    r._parents = nps
    r.parents = set([obodag[p] for p in nps])

for go, r in obodag.items():
    if go in bad_go:
        del(obodag[go])


for g in geneid2gos:
    geneid2gos[g] = set(go for go in geneid2gos[g] if go in obodag)


reverse = {}
for k in GeneID2nt:
    reverse[GeneID2nt[k].Symbol] = k


if 'apv' in strategy:
    which = 5  # 5 is adj pv, 4 is pv
else:
    which = 4  # 5 is adj pv, 4 is pv

if '0.' in strategy:
    pv_cut = float('0.' + strategy.split('0.')[1])

GO_results = {}

if strategy == 'switch':
    test1 = lambda s, x, y: results[s][x,y][0] <= 0.5 and results[s][x,y][1] >= 0.5
    test2 = lambda s, x, y: results[s][x,y][0] >= 0.5 and results[s][x,y][1] <= 0.5
elif 'ipv' in  strategy:
    test1 = lambda s, x, y: results[s][x,y][which] > pv_cut and results[s][x,y][0] > .7 and results[s][x,y][1] > .7
    test2 = lambda s, x, y: results[s][x,y][which] > pv_cut and results[s][x,y][0] < .3 and results[s][x,y][1] < .3
elif 'Mpv' in  strategy:
    test1 = lambda s, x, y: results[s][x,y][which] <= pv_cut and results[s][x,y][3] > 0 and results[s][x,y][0] > .5
    test2 = lambda s, x, y: results[s][x,y][which] <= pv_cut and results[s][x,y][3] < 0 and results[s][x,y][0] < .5
elif 'Lpv' in  strategy:
    test1 = lambda s, x, y: results[s][x,y][which] <= pv_cut and results[s][x,y][3] > 0 and results[s][x,y][1] < .5
    test2 = lambda s, x, y: results[s][x,y][which] <= pv_cut and results[s][x,y][3] < 0 and results[s][x,y][1] > .5
elif 'pv' in  strategy:
    test1 = lambda s, x, y: results[s][x,y][which] <= pv_cut and results[s][x,y][3] > 0
    test2 = lambda s, x, y: results[s][x,y][which] <= pv_cut and results[s][x,y][3] < 0

for spe in species:
    goeaobj = GOEnrichmentStudy(
            [k for k in GeneID2nt],  # we take ALL GENES
            geneid2gos, # geneid/GO associations
            obodag, # Ontologies
            propagate_counts = False,
            verbose=False,
            alpha = 0.05, # default significance cut-off
            methods = ['fdr_bh']) # defult multipletest correction method
    if from_anc == 'primate' and spe != 'Mus musculus':
        key = (spe, '-'.join(sorted((tree & ('Callithrix jacchus' if from_anc == 'primate' else spe)
                                    ).up.get_leaf_names())))
    elif from_anc == 'anc':
        key = (spe, '-'.join(sorted((tree & (spe)).up.get_leaf_names())))
    else:
        if spe == 'Homo sapiens':
            continue
        key = ('Homo sapiens', spe)
    GO_results[spe] = {}
    print (spe)
    print ('.' * 50)

    pos_list = [(c, b) for c, b in results[key]
                if test1(key, c, b) and np.isfinite(results[key][c, b][which])]
    gene_list = reduce(lambda x, y: x + y, [genes.get((c, b), []) for c, b in pos_list])

    # write gene_list
    out = open('coord_{}_A-{}_B-{}.tsv'.format(strategy, species[key[0]]['vulgar'], species[key[1]]['vulgar']), 'w')
    out.write('\n'.join(['{}\t{}\t{}'.format(c, b*reso, b * reso + reso) for c, b in pos_list]) + '\n')
    out.close()
    out = open('genes_{}_A-{}_B-{}.tsv'.format(strategy, species[key[0]]['vulgar'], species[key[1]]['vulgar']), 'w')
    out.write('\n'.join(gene_list) + '\n')
    out.close()

    # run GSE test
    print (' - Repressed ({:5,d}): '.format(len(gene_list)))
    results_all = goeaobj.run_study([reverse[k] for k in gene_list if k in reverse], log=None)
    results_sig = [r for r in results_all if r.p_fdr_bh < 0.05]
    GO_results[spe]['Repressed'] = results_all
    if results_sig:
        for r, s in sorted([(r.get_pvalue(), r.name) for r in results_sig]):
            print ('   -> {:8.3g} {}'.format(r, s))
        plot_results("sign_{}_pv05_repressed.png".format(spe.replace(' ', ' ')), results_sig,
                     log=open('/dev/null', 'w'))
    print ('.' * 50)

    pos_list = [(c, b) for c, b in results[key]
                if test2(key, c, b) and np.isfinite(results[key][c, b][which])]
    gene_list = reduce(lambda x, y: x + y, [genes.get((c, b), []) for c, b in pos_list])

    # write gene_list
    out = open('coord_{}_B-{}_A-{}.tsv'.format(strategy, species[key[0]]['vulgar'], species[key[1]]['vulgar']), 'w')
    out.write('\n'.join(['{}\t{}\t{}'.format(c, b*reso, b * reso + reso) for c, b in pos_list]) + '\n')
    out.close()
    out = open('genes_{}_B-{}_A-{}.tsv'.format(strategy, species[key[0]]['vulgar'], species[key[1]]['vulgar']), 'w')
    out.write('\n'.join(gene_list) + '\n')
    out.close()

    # run GSE test
    print (' - Activated ({:5,d}): '.format(len(gene_list)))
    results_all = goeaobj.run_study([reverse[k] for k in gene_list if k in reverse], log=None)
    results_sig = [r for r in results_all if r.p_fdr_bh < 0.05]
    GO_results[spe]['Activated'] = results_all
    if results_sig:
        for r, s in sorted([(r.get_pvalue(), r.name) for r in results_sig]):
            print ('   -> {:8.3g} {}'.format(r, s))
        plot_results("sign_{}_pv05_activated.png".format(spe.replace(' ', ' ')), results_sig,
                     log=open('/dev/null', 'w'))
    print ('=' * 50)


dots = {'Repressed': {}, 'Activated': {}}
for nspe, spe in enumerate(sorted(GO_results.keys(), key=lambda x: species.keys().index(x))):
    for state in dots:
        for r in GO_results[spe][state]:
            if r.p_fdr_bh > 0.05:
                continue
            dots[state].setdefault(r.name, [None] * len(species))
            dots[state][r.name][nspe] = r


all_names = list(set(dots['Activated'].keys() + dots['Repressed'].keys()))


go2name = dict([(r.GO, n) for r in dots['Activated'].get(n, dots['Repressed'].get(n, '')) if r][0]
               for n in all_names)
name2go = dict([(n, r.GO) for r in dots['Activated'].get(n, dots['Repressed'].get(n, '')) if r][0]
               for n in all_names)

termcount = TermCounts(obodag, geneid2gos)


dist_matrix = [[((lin_sim(name2go[go1], name2go[go2], obodag, termcount) if go1 != go2 else 0) or 0)
                for go1 in all_names] for go2 in all_names]
dist_cond = [dist_matrix[i][j] for i in range(len(dist_matrix)) for j in range(i + 1, len(dist_matrix))]


dist_matrix = [[(1 + dist_matrix[i][j]) if i!= j else 0
                for i in range(len(dist_matrix))] for j in range(len(dist_matrix))]


for i in range(len(dist_matrix)):
    for j in range(i + 1, len(dist_matrix)):
        if dist_matrix[i][j] != dist_matrix[j][i]:
            print i, j, dist_matrix[i][j], dist_matrix[j][i]
            dist_matrix[i][j] = dist_matrix[j][i] = min((dist_matrix[i][j], dist_matrix[j][i]))

dist_cond = squareform(dist_matrix)

plt.figure(figsize=(18, 15))
plt.imshow(dist_matrix, cmap='viridis_r')
plt.colorbar()
_= plt.yticks(range(len(all_names)), all_names)


link = linkage(dist_cond, method='ward')


plt.figure(figsize=(4, 17))
dd = dendrogram(link, labels=all_names, orientation='right', leaf_font_size=12)#, link_color_func=lambda x: 'grey')


minmax = lambda x: (min(x), max(x))


sorted_names = dd['ivl']


idx = [all_names.index(n) for n in sorted_names]


dist_matrix = [[dist_matrix[i][j] for j in idx] for i in idx]


colors = ['#ff7f0e', '#9467bd', '#2ca02c', '#1f77b4', '#d62728', '#8c564b',
          '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']

cl_coords1 = [i[0] for i in cut_tree(link, n_clusters=None, height=0.9)[idx]]
cl_coords2 = [i[0] for i in cut_tree(link, n_clusters=3)[idx]]
cl_coords, cl_colors = zip(*[(c1, colors[c2]) for c1, c2 in zip(cl_coords1, cl_coords2)])
cl_colors = [cl_colors[[n for n, c in enumerate(cl_coords) if c == i][0]] for i in range(len(set(cl_coords)))]
cl_coords = [minmax([n for n, c in enumerate(cl_coords) if c == i]) for i in range(len(set(cl_coords)))]


# In[ ]:


f = plt.figure(figsize=(18, 15))
axe = f.add_subplot(111)
plt.imshow(dist_matrix, cmap='viridis_r')
plt.colorbar()
_= plt.yticks(range(len(sorted_names)), [n + '   ' for n in sorted_names])
for i, (y1, y2) in enumerate(cl_coords):
    axe.add_patch(plt.Rectangle((-1.5, y1-0.45), 0.8, y2-y1 + 0.9, clip_on=False,
                                facecolor=cl_colors[i], alpha=0.6))


rpath = '/data/Projects/GEVO_3D/results/compartments/functional_enrichment/'


fig = plt.figure(figsize=(12, 0.18 * len(all_names) + 3))

x, y, c, z, R = zip(*[(x, y,
                       r.p_fdr_bh if r else 1,
                       r.ratio_in_study[0] if r else 0,
                       (np.divide(*map(float, r.ratio_in_study)) / np.divide(*map(float, r.ratio_in_pop))
                        if r else 0))
                    for y, n in enumerate(sorted_names)
                    for x, r in enumerate(dots['Activated'].get(n, [0,0,0,0,0,0,0,0]))])

z = [np.log(v) * 50 if v else 0 if c[i] == 1 else 25 for i, v in enumerate(z)]
c = [(1 if r > 1 else -1) * v for v, r in zip(c, R)]

axe1 = fig.add_axes([0.56, 0.2, 0.15, 0.75])
if 'ipv' in strategy:
    axe1.set_title('Most conserved\nA compartments', size=11)
else:
    axe1.set_title('Less active in {}'.format('human' if from_anc == 'nop' else 'ancestor'), size=11)
sc1 = axe1.scatter(x, y, s=z, c=[(1 if v > 0 else -1) * np.log(1/abs(v)) for v in c],
                   alpha=0.8, cmap='coolwarm', vmin=-10, vmax=10)
axe1.grid(alpha=0.3)
_ = axe1.set_yticks(range(len(sorted_names)))
_ = axe1.set_yticklabels([n + ' ' for n in sorted_names])
_ = axe1.set_xticks(range(len(GO_results)))
_ = axe1.set_xticklabels([species[spe]['vulgar'] for spe in species
                         if from_anc != 'nop' or spe != 'Homo sapiens'], rotation=-90,
                         size=13, ha='center')
axe1.set_axisbelow(True)
axe1.set_xlim(-0.75, len(GO_results) - .25)

x, y, c, z, R = zip(*[(x, y,
                       r.p_fdr_bh if r else 1,
                       r.ratio_in_study[0] if r else 0,
                       (np.divide(*map(float, r.ratio_in_study)) / np.divide(*map(float, r.ratio_in_pop))
                        if r else 0))
                      for y, n in enumerate(sorted_names)
                      for x, r in enumerate(dots['Repressed'].get(n, [0,0,0,0,0,0,0]))])

z = [np.log(v) * 50 if v else 0 if c[i] == 1 else 25 for i, v in enumerate(z)]
c = [(1 if r > 1 else -1) * v for v, r in zip(c, R)]

axe2 = fig.add_axes([0.72, 0.2, 0.15, 0.75], sharey=axe1)
if 'ipv' in strategy:
    axe2.set_title('Most conserved\nB compartments', size=11)
else:
    axe2.set_title('More active in {}'.format('human' if from_anc == 'nop' else 'ancestor'), size=11)
sc2 = axe2.scatter(x, y, s=z, c=[(1 if v > 0 else -1) * np.log(1/abs(v)) for v in c],
                   alpha=0.8, cmap='coolwarm', vmin=-10, vmax=10)
axe2.grid(alpha=0.3)
# _ = axe2.set_yticks(range(len(sorted_names)))
# _ = axe2.set_yticklabels([])
plt.setp(axe2.get_yticklabels(), visible=False)
_ = axe2.set_xticks(range(len(species)))
_ = axe2.set_xticklabels([species[spe]['vulgar'] for spe in species
                          if from_anc != 'nop' or spe != 'Homo sapiens'], rotation=-90,
                         size=13, ha='center')
axe2.yaxis.set_ticks_position('none')
axe2.set_axisbelow(True)
axe2.set_xlim(-0.75, len(GO_results) - .25)
axe2.set_ylim((-0.9, len(sorted_names) - 0.1))

cax = fig.add_axes([0.885, 0.25,0.02,0.3])
cb = plt.colorbar(sc1, cax=cax)
ticks = [-0.1, -0.01, -0.001, -0.0001, 0, 0.0001, 0.001, 0.01, 0.1]
cb.set_ticks([(1 if t > 0 else -1) * np.log(1. / abs(t)) if t else 0 for t in ticks])
cb.set_ticklabels([abs(t) for t in ticks])
cax.set_ylabel('Adjusted p-value\nOver-represented / Under-represented', rotation=-90, va='bottom')

g0 = plt.scatter([],[], s=25               , marker='o', color='', edgecolor='k')
g1 = plt.scatter([],[], s=np.log(10  ) * 50, marker='o', color='', edgecolor='k')
g2 = plt.scatter([],[], s=np.log(100 ) * 50, marker='o', color='', edgecolor='k')
g3 = plt.scatter([],[], s=np.log(1000) * 50, marker='o', color='', edgecolor='k')

_ = axe2.legend((g0, g1, g2, g3), ('0', '10', '100', '1000'),
                title='Number of genes', frameon=False,
                scatterpoints=1, labelspacing=0.7,
                bbox_to_anchor=(1, 1), ncol=1)

for i, (y1, y2) in enumerate(cl_coords):
    axe1.add_patch(plt.Rectangle((-1, y1-0.45), 0.3, y2-y1 + 0.9, clip_on=False,
                                 facecolor=cl_colors[i], alpha=0.6))
plt.savefig(os.path.join(rpath, 'Goterms_cutlvl%d_%s_%s_from-%s.pdf' % (cutlvl, strategy, goset, from_anc)),
            format='pdf')
