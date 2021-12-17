import os
import re
from argparse     import ArgumentParser


def _find_region(query, beg_pos, end_pos, sam_handler, pattern):
    new_pos = int((end_pos + beg_pos) / 2 )
    sam_handler.seek(new_pos)
    next(sam_handler)
    line = next(sam_handler)
    new_read = int(pattern.findall(line)[0])
    # print(query, new_read, beg_pos, end_pos)
    if new_read > query:    # almost there, try a bigger number
        # print('try bigger')
        end_pos = new_pos
    elif new_read < query:  # almost there, try a smaller number
        # print('try smaller')
        beg_pos = new_pos
    else:                   # congratulations you won!!!!
        return line
    return _find_region(query, beg_pos, end_pos, sam_handler, pattern)


def find_region(query, sam_file, pattern):
    query = int(pattern.findall(query)[0])
    sam_handler = open(sam_file)

    # get start
    beg_pos = 0
    for line in sam_handler:
        if not line.startswith('@'):
            beg_pos -= len_line
            break
        len_line = len(line)
        beg_pos += len_line

    # get end
    end_pos = os.path.getsize(sam_file)

    return(_find_region(query, beg_pos, end_pos, sam_handler, pattern))

def main():
    opts = get_options()
    pattern = re.compile(opts.pattern)
    
    region = find_region(opts.query, opts.sam_file, pattern)

    print(region, end='')


def get_options():
    parser = ArgumentParser()

    parser.add_argument('-i', dest='sam_file', metavar='PATH', required=True,
                        default=False, help='input SAM file sorted by QNAME.')
    parser.add_argument('-q', '--query', dest='query', required=True,
                        help='''Query name from which we want to extract data 
                        from the SAM.''')
    parser.add_argument('--pattern', dest='pattern', default="^([0-9]+)_",
                        help='''[%(default)s] regexp pattern to extract number 
                        (on which sorting is made) from QNAME.''')
    opts = parser.parse_args()
    return opts

if __name__ == "__main__":
    exit(main())