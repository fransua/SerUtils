"""
"""

from matplotlib         import pyplot as plt
from matplotlib.patches import FancyBboxPatch
from matplotlib.ticker  import FuncFormatter

def read_genes(fnam):
    """
    Parse file from biomart-ensembl

    downloaded using this XML:

    xml = '''
    <?xml version="1.0" encoding="UTF-8"?>
    <!DOCTYPE Query>
    <Query  virtualSchemaName = "default" formatter = "TSV" header = "1" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" >

    	<Dataset name = "hsapiens_gene_ensembl" interface = "default" >
    		<Attribute name = "ensembl_gene_id" />
    		<Attribute name = "ensembl_transcript_id" />
    		<Attribute name = "strand" />
    		<Attribute name = "chromosome_name" />
    		<Attribute name = "start_position" />
    		<Attribute name = "end_position" />
    		<Attribute name = "cds_start" />
    		<Attribute name = "cds_end" />
    		<Attribute name = "exon_chrom_start" />
    		<Attribute name = "exon_chrom_end" />
    		<Attribute name = "external_gene_name" />
    		<Attribute name = "gene_biotype" />
    	</Dataset>
    </Query>
    '''
    xml = ''.join(l.strip() for l in xml.split('\n'))

    ! wget -O result.txt 'http://www.ensembl.org/biomart/martservice?query={xml}'

    :returns gene disctionary containting longest transcripts only:
    """
    fh = open(fnam, 'r')
    header = next(fh)
    header = header.strip().split('\t')
    vtype = {'Gene stable ID': str,
             'Transcript stable ID': str,
             'Strand': str,
             'Chromosome/scaffold name': str,
             'Gene start (bp)': int,
             'Gene end (bp)': int,
             'CDS start': int,
             'CDS end': int,
             'Exon region start (bp)': int,
             'Exon region end (bp)': int,
             'Gene name': str,
             'Gene type': str}
    transcripts = {}
    for line in fh:
        tmp = dict(zip(header, line.strip().split('\t')))
        for k, v in tmp.items():
            try:
                tmp[k] = vtype[k](v)
            except ValueError:
                tmp[k] = float('nan')
        transcripts.setdefault(tmp['Transcript stable ID'], []).append(tmp)
    genes = {}
    for t in transcripts:
        val = sum(tmp['Exon region end (bp)'] - tmp['Exon region start (bp)'] for tmp in transcripts[t])
        if val <= 0:
            break
        genes.setdefault(transcripts[t][0]['Gene stable ID'], []).append((t, val))
    # get longest transcript
    for g in genes:
        t = max(genes[g], key=lambda x: x[1])[0]
        genes[g] = t

    good_transcripts = set(genes.values())
    del(genes)

    transcripts = dict((t, transcripts[t]) for t in transcripts if t in good_transcripts)

    for t in good_transcripts:
        exons = [(tmp['Exon region start (bp)'], tmp['Exon region end (bp)']) for tmp in transcripts[t]]
        transcripts[t] = dict((k, transcripts[t][0][k])
                              for k in ['Gene stable ID', 'Strand', 'Chromosome/scaffold name',
                                        'Gene start (bp)', 'Gene end (bp)', 'CDS start', 'CDS end',
                                        'Gene name', 'Gene type'])
        transcripts[t]['exons'] = sorted(exons)
        # some filtering...
        # if len(exons) == 1 and (exons[0][1] - exons[0][0]) < 333:
        #     del transcripts[t]
    return transcripts


def nicer(res, sep=' ', comma='', allowed_decimals=0):
    """
    writes resolution number for human beings.

    :param ' ' sep: character between number and unit (e.g. default: '125 kb')
    :param '' comma: character to separate groups of thousands
    :param 0 allowed_decimals: if 1 '1900 kb' would be written as '1.9 Mb'
    """
    format = lambda x: '{:,g}'.format(x).replace(',', comma)

    if not res:
        return format(res) + sep + 'b'
    if not res % 10**(9 - allowed_decimals):
        return format(res / 10.**9) + sep + 'Gb'
    if not res % 10**(6 - allowed_decimals):
        return format(res / 10.**6) + sep + 'Mb'
    if not res % 10**(3 - allowed_decimals):
        return format(res / 10.**3) + sep + 'kb'
    return format(res) + sep + 'b'


def get_text_width(txt_obj, renderer, axe):
    bb = txt_obj.get_window_extent(renderer=renderer)

    x, y = axe.transData.transform_point((3, 4))

    inv = axe.transData.inverted()
    (x1, y1), (x2, y2) = inv.transform(((x, y), (x + bb.width, y + bb.height)))
    return x2 - x1


def plot_genes(genes, crm, beg, end, gene_type_color=None, axe=None,
               extra_fig_heigth=0, highlight=None):
    """
    plot genes through their exons/introns.

    :param genes: dictionary returned by read_genes function
    :param crm: Chromosome name of the wanted region
    :param beg: start position of the wanted region
    :param end: end position of the wanted region
    :param None axe: Matplotlib Axe object
    """
    if not gene_type_color:
        gene_type_color = {'protein_coding'                     : 'firebrick',
                           'polymorphic'                        : 'darkred',
                           'processed_transcript'               : 'forestgreen',
                           'lncRNA'                             : 'forestgreen',
                           'non_coding'                         : 'forestgreen',
                           '3prime_overlapping_ncRNA'           : 'forestgreen',
                           'antisense'                          : 'forestgreen',
                           'lincRNA'                            : 'forestgreen',
                           'retained_intron'                    : 'forestgreen',
                           'sense_intronic'                     : 'forestgreen',
                           'sense_overlapping'                  : 'forestgreen',
                           'macro_lncRNA'                       : 'forestgreen',
                           'ncRNA'                              : 'green',
                           'miRNA'                              : 'green',
                           'misc_RNA'                           : 'darkgreen',
                           'piRNA'                              : 'green',
                           'rRNA'                               : 'green',
                           'sRNA'                               : 'green',
                           'scRNA'                              : 'green',
                           'scaRNA'                             : 'green',
                           'siRNA'                              : 'green',
                           'snRNA'                              : 'green',
                           'snoRNA'                             : 'green',
                           'tRNA'                               : 'green',
                           'ribozyme'                           : 'green',
                           'Mt_rRNA'                            : 'green',
                           'Mt_tRNA'                            : 'green',
                           'vaultRNA'                           : 'darkgreen',
                           'rRNA_pseudogene'                    : 'darkgreen',
                           'unclassified_processed_transcript'  : 'steelblue',
                           'pseudogene'                         : 'steelblue',
                           'processed_pseudogene'               : 'steelblue',
                           'unprocessed_pseudogene'             : 'steelblue',
                           'bidirectional_promoter_lncRNA'      : 'steelblue',
                           'transcribed_pseudogene'             : 'steelblue',
                           'transcribed_unprocessed_pseudogene' : 'steelblue',
                           'transcribed_processed_pseudogene'   : 'steelblue',
                           'transcribed_unitary_pseudogene'     : 'steelblue',
                           'translated_pseudogene'              : 'steelblue',
                           'translated_processed_pseudogene'    : 'steelblue',
                           'polymorphic_pseudogene'             : 'steelblue',
                           'unitary_pseudogene'                 : 'steelblue',
                           'IG_pseudogene'                      : 'steelblue',
                           'IG_C_pseudogene'                    : 'steelblue',
                           'IG_V_pseudogene'                    : 'steelblue',
                           'IG_J_pseudogene'                    : 'steelblue',
                           'IG_gene'                            : 'darkred',
                           'IG_C_gene'                          : 'darkred',
                           'IG_D_gene'                          : 'darkred',
                           'IG_V_gene'                          : 'darkred',
                           'IG_J_gene'                          : 'darkred',
                           'TR_pseudogene'                      : 'steelblue',
                           'TR_C_pseudogene'                    : 'steelblue',
                           'TR_D_pseudogene'                    : 'steelblue',
                           'TR_V_pseudogene'                    : 'steelblue',
                           'TR_J_pseudogene'                    : 'steelblue',
                           'TR_gene'                            : 'darkred',
                           'TR_C_gene'                          : 'darkred',
                           'TR_D_gene'                          : 'darkred',
                           'TR_V_gene'                          : 'darkred',
                           'TR_J_gene'                          : 'darkred',
                           'TEC'                                : 'grey'}
    if not axe:
        fig = plt.figure(figsize=(10, 3))
        axe = plt.subplot(111)
        renderer = fig.canvas.get_renderer()
    else:
        fig = axe.get_figure()
        renderer = fig.canvas.get_renderer()

    highlight = set(highlight) if highlight else set()
    axe.set_xlim(beg, end)

    minimum = float('inf')
    maximum = float('-inf')
    prev_g, prev_e, prev_s = None, None, None

    nrows = 1

    min_dist = (end - beg) / 500

    pos = []
    for g in sorted([g for g in genes if
                     genes[g]['Chromosome/scaffold name'] == crm and
                     genes[g]['exons'][0][0] >= beg and
                     genes[g]['exons'][-1][1] <= end],
                    key=lambda x: genes[x]['exons'][0][0]):

        begs = tuple([b for b, _ in genes[g]['exons']])
        vs = sorted([b for b, _ in genes[g]['exons']] + [e for _, e in genes[g]['exons']])
        b = vs[0]
        e = vs[-1]
        seen = []
        for bb, ee, nnum in pos:
            if bb <= b <= ee or bb <= e <= ee:
                seen.append(nnum)
            elif e < bb:
                break
        try:
            num = sorted(list(set(seen) ^ set(range(0, max(seen) + 2))))[0]
        except ValueError:
            num = 0
        color = gene_type_color.get(genes[g]['Gene type'], 'red')
        if g in highlight:
            t = axe.text(genes[g]['exons'][0][0],
                     num + 0.1,
                     genes[g]['Gene name'] + ('>' if genes[g]['Strand'] == '1' else '<'),
                         color=color, size=9, weight="bold")
        else:
            t = axe.text(genes[g]['exons'][0][0],
                     num + 0.1,
                     genes[g]['Gene name'] + ('>' if genes[g]['Strand'] == '1' else '<'),
                         color=color, size=8)
        w = get_text_width(t, renderer, axe)
        pos.append((b, max(b + w, e), num))
        for b, e in genes[g]['exons']:
            if prev_g == g:
                if prev_e and min_dist < b - prev_e:
                    axe.plot([prev_e, (b + prev_e) / 2., b],
                             [num - 0.2, num - 0.05, num - 0.2], color=color, lw=1)
            p_fancy = FancyBboxPatch((b, num - 0.25), e - b, 0.15, alpha=1,
                                     color=color)
            axe.add_patch(p_fancy)
            p_fancy.set_boxstyle("Round", pad=0.1, rounding_size=0.3)

            prev_g, prev_e = g, e
            minimum = min(b, minimum)
            maximum = max(e, maximum)
        try:
            nrows = max(num + 1, nrows)
        except ValueError:
            pass
    axe.set_ylim(-0.5, nrows - 0.5)
    axe.spines['right'].set_visible(False)
    axe.spines['left'].set_visible(False)
    axe.spines['top'].set_visible(False)
    axe.set_yticks([])
    axe.xaxis.grid(alpha=0.3, lw=1.5)

    def format_xticks(tickstring, _=None):
        tickstring = int(tickstring + beg)
        return nicer(tickstring if tickstring else 1,
                     comma=',', allowed_decimals=1, sep='')

    axe.xaxis.set_major_formatter(FuncFormatter(format_xticks))

    labels = axe.get_xticklabels()
    plt.setp(labels, rotation=0, ha='center', size=10)

    xlabel = '{}: {:,}-{:,}'.format(('' if 'chr' in crm else 'chr') + crm, beg if beg else 1, end)
    axe.set_xlabel(xlabel)
    fig.set_figheight(0.3 + nrows / 3. + extra_fig_heigth)
    return nrows
