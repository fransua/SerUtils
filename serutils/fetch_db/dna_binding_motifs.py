"""
"""
import xml.etree.ElementTree as ET
from re import findall


def reverse_factorbook_seq_logo(fname):
    """
    from an SVG file from https://factorbook.org/

    :returns: a matrix woth proportions
    """
    tree = ET.parse(fname)
    root = tree.getroot()
    matrix = []
    for elem in root.iter():
        data = elem.get('transform')
        if not data:
            continue
        if (data.startswith('translate(') and not data.startswith('translate(0,') and
            not data.endswith('(0,0)')):
            matrix.append({})
            continue
        elif data.startswith('scale'):
            cur_scale = float(findall('scale\([0-9.]+,([0-9.]+)\)', data)[0])
            letter = elem.getchildren()[0].get('regex')
            if not letter:
                letter = elem.getchildren()[0].getchildren()[0].get('regex')
            matrix[-1][letter] = cur_scale

    todel = []
    for i, col in enumerate(matrix):
        if not col:
            todel.append(i)
            continue
        tot = sum(col.values())
        for let in 'ACGT':
            try:
                col[let] /= tot
            except KeyError:
                col[let] = 0.0
    for i in todel[::-1]:
        del(matrix[i])
    return matrix
