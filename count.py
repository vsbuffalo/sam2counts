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
    """
    counts = dict()
    sf = pysam.Samfile(filename, 'r')
    for read in sf:
        id_name = sf.getrname(read.rname) if read.rname != -1 else 0
        if not read.is_unmapped:
            counts[id_name] = counts.get(id_name, 0) + 1
        elif read.is_unmapped:
            # we want to record these as zeros
            counts[id_name] = counts.get(id_name, 0)
    return counts

def collapsed_nested_count_dict(counts_dict, all_ids, order=None):
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
    if order is None:
        col_order = total_counts.keys()

    collapsed_dict = dict()
    for i, filename in enumerate(total_counts.keys()):
        for id_name in all_ids:
            if not collapsed_dict.get(id_name, False):
                collapsed_dict[id_name] = list()

            # get counts and append
            c = total_counts[filename].get(id_name, 0)
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

    parser.add_option("-v", "--verbose", dest="verbose",
                      help="enable verbose output")

    (options, args) = parser.parse_args()

    if len(args) < 1:
        parser.error("one or more SAM files as arguments required")

    total_counts = dict()
    all_ids = list()
    for filename in args:
        total_counts[filename] = SAM_file_to_counts(filename)
        all_ids.extend(total_counts[filename].keys())

    all_ids = set(all_ids)
    table_dict = collapsed_nested_count_dict(total_counts, all_ids)
    counts_to_file(table_dict, 'test.csv', delimiter=options.delimiter)
        
        
    
    
    
