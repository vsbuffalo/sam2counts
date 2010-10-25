"""
count.py -- Take SAM files and output a table of counts with column
names that are the filenames, and rowname that are the reference
names.

Author: Vince Buffalo
Email: vsbuffaloAAAAA@gmail.com (with poly-A tail removed)
"""

import sys
import csv
try:
    import pysam
except ImportError:
    sys.exit("pysam not installed; please install it\n")

from optparse import OptionParser

def SAM_file_to_counts(filename):
    """
    Take a filename to a SAM file, and create a hash of mapped and
    unmapped reads; values are the counts of occurences.

    Also, a hash of qualities (either 0 or otherwise) of mapped reads
    is output.
    """
    counts = dict()
    qual_counts = dict()
    sf = pysam.Samfile(filename, 'r')
    for read in sf:
        id_name = sf.getrname(read.rname) if read.rname != -1 else 0

        ## quality recording
        if not qual_counts.get(id_name, False):
            # create empty list for hash
            qual_counts[id_name] = {'0':0, '>0':0}

        ## initiate entry; even if not mapped, record 0 count
        counts[id_name] = counts.get(id_name, 0)
        
        if not read.is_unmapped:
            counts[id_name] = counts.get(id_name, 0) + 1

            if read.mapq == 0:
                qual_counts[id_name]['0'] = qual_counts[id_name]['0'] + 1
            else:
                qual_counts[id_name]['>0'] = qual_counts[id_name]['>0'] + 1

    return {'counts':counts, 'qual_counts':qual_counts}

def collapsed_nested_count_dict(counts_dict, all_ids):
    """
    This function takes a nested dictionary `counts_dict` and
    `all_ids`, which is built with the `table_dict`. All files (first
    keys) in `counts_dict` are made into columns with order specified
    by `order`.

    Output is a dictionary with keys that are the id's (genes or
    transcripts), with values that are ordered counts. A header will
    be created on the first row from the ordered columns (extracted
    from filenames).
    """
    col_order = counts_dict.keys()

    collapsed_dict = dict()
    for i, filename in enumerate(counts_dict.keys()):
        for id_name in all_ids:
            if not collapsed_dict.get(id_name, False):
                collapsed_dict[id_name] = list()

            # get counts and append
            c = counts_dict[filename].get(id_name, 0)
            collapsed_dict[id_name].append(c)
    return {'table':collapsed_dict, 'header':col_order}


def counts_to_file(table_dict, outfilename, delimiter=','):
    """
    A function for its side-effect of writing `table_dict` (which
    contains a table and header), to `outfilename` with the specified
    `delimiter`.
    """
    writer = csv.writer(open(outfilename, 'w'), delimiter=delimiter)
    table = table_dict['table']
    header = table_dict['header']
    header_row = True
    
    for id_name, fields in table.items():
        if header_row:
            row = ['id'] + header
            writer.writerow(row)
            header_row = False

        if id_name == 0:
            continue
        row = [id_name]
        row.extend(fields)
        writer.writerow(row)

if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option("-d", "--delimiter", dest="delimiter",
                      help="the delimiter (default: tab)", default='\t')
    parser.add_option("-o", "--out-file", dest="out_file",
                      help="the output file's name (default: counts.txt)",
                      default='counts.txt')
    parser.add_option("-q", "--include-quality", dest="quality",
                      help="generate quality statistics too (default: False)",
                      default=False)
    parser.add_option("-v", "--verbose", dest="verbose",
                      help="enable verbose output")

    (options, args) = parser.parse_args()

    if len(args) < 1:
        parser.error("one or more SAM files as arguments required")

    file_counts = dict()
    file_qual_counts = dict()
    all_ids = list()
    for filename in args:
        ## read in SAM file, extract counts, and unpack counts and qual_counts
        tmp = SAM_file_to_counts(filename)
        counts, qual_counts = tmp['counts'], tmp['qual_counts']

        ## save counts, qual_counts, and all ids encountered
        file_counts[filename] = counts
        file_qual_counts[filename] = qual_counts
        all_ids.extend(file_counts[filename].keys())

    ## Uniquify all_ids, and then take the nested file_counts
    ## dictionary, collapse, and write to file.
    all_ids = set(all_ids)
    table_dict = collapsed_nested_count_dict(file_counts, all_ids)
    counts_to_file(table_dict, options.out_file, delimiter=options.delimiter)
        
        
    
    
    
