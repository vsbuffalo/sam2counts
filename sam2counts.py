"""
count.py -- Take SAM files and output a table of counts with column
names that are the filenames, and rowname that are the reference
names.

Author: Vince Buffalo
Email: vsbuffaloAAAAA@gmail.com (with poly-A tail removed)
"""

VERSION = 0.91

import sys
import csv
from os import path
try:
    import pysam
except ImportError:
    sys.exit("pysam not installed; please install it\n")

from optparse import OptionParser

# csv module stupidly uses \r\n by default; this is just a Unix
# variant of the "Excel" dialect.
UnixVariant = csv.excel
UnixVariant.lineterminator = "\n"

def SAM_file_to_counts(filename, bam=False, extra=False, use_all_references=True, min_mapq=None):
    """
    Take SAM filename, and create a hash of mapped and unmapped reads;
    keys are reference sequences, values are the counts of occurences.

    Also, a hash of qualities (either 0 or >0) of mapped reads
    is output, which is handy for diagnostics.
    """
    counts = dict()
    unique = dict()
    nonunique = dict()
    mode = 'r'
    if bam:
        mode = 'rb'
    sf = pysam.Samfile(filename, mode)

    if use_all_references:
        # Make dictionary of all entries in header
        for sn in sf.header['SQ']:
            if extra:
                unique[sn['SN']] = 0
                nonunique[sn['SN']] = 0
            counts[sn['SN']] = 0

    for read in sf:
        if not read.is_unmapped:
            if min_mapq is not None and read.mapq < min_mapq:
                continue
            id_name = sf.getrname(read.rname) if read.rname != -1 else 0

            if not use_all_references and not counts.get(id_name, False):
                ## Only make keys based on aligning reads, make empty hash
                if extra:
                    unique[id_name] = 0
                    nonunique[id_name] = 0
                ## initiate entry; even if not mapped, record 0 count
                counts[id_name] = counts.get(id_name, 0)
        
            
            counts[id_name] = counts.get(id_name, 0) + 1

            if extra:
                if read.mapq == 0:
                    nonunique[id_name] = nonunique[id_name] + 1
                else:
                    unique[id_name] = unique[id_name] + 1

    if extra:
        return {'counts':counts, 'unique':unique, 'nonunique':nonunique}

    return {'counts':counts}

def collapsed_nested_count_dict(counts_dict, all_ids, order=None):
    """
    Takes a nested dictionary `counts_dict` and `all_ids`, which is
    built with the `table_dict`. All files (first keys) in
    `counts_dict` are made into columns with order specified by
    `order`.

    Output is a dictionary with keys that are the id's (genes or
    transcripts), with values that are ordered counts. A header will
    be created on the first row from the ordered columns (extracted
    from filenames).
    """
    if order is None:
        col_order = counts_dict.keys()
    else:
        col_order = order

    collapsed_dict = dict()
    for i, filename in enumerate(col_order):
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
    writer = csv.writer(open(outfilename, 'w'), delimiter=delimiter, dialect=UnixVariant)
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
    parser.add_option("-m", "--min-mapq", type=int,
                      help="minimum mapping quality (default: none)", default=None)
    parser.add_option("-o", "--out-file", dest="out_file",
                      help="output filename (default: counts.txt)",
                      default='counts.txt', action="store", type="string")
    parser.add_option("-u", "--extra-output", dest="extra_out",
                      help="output extra information on non-unique and unique mappers (default: False)",
                      default=False, action="store_true")
    parser.add_option("-b", "--bam", dest="bam",
                      help="all input files are BAM (default: False)", 
                      default=False, action="store_true")
    parser.add_option("-r", "--use-all-references", dest="use_all_references",
                      help="Use all the references from the SAM header (default: True)",
                      default=True, action="store_false")
    parser.add_option("-f", "--extra-out-files", dest="extra_out_files",
                      help="comma-delimited filenames of unique and non-unique output "
                      "(default: unique.txt,nonunique.txt)",
                      default='unique.txt,nonunique.txt', action="store", type="string")
    parser.add_option("-v", "--verbose", dest="verbose",
                      help="enable verbose output")

    (options, args) = parser.parse_args()

    if len(args) < 1:
        parser.error("one or more SAM files as arguments required")

    file_counts = dict()
    file_unique_counts = dict()
    file_nonunique_counts = dict()
    all_ids = list()
    files = [path.basename(f) for f in args]

    if len(set(files)) != len(set(args)):
        parser.error("file args must have unique base names (i.e. no foo/bar joo/bar)")

    ## do a pre-run check that all files exist
    for full_filename in args:
        if not path.exists(full_filename):
            parser.error("file '%s' does not exist" % full_filename)
        
    for full_filename in args:
        filename = path.basename(full_filename)
        ## read in SAM file, extract counts, and unpack counts
        tmp = SAM_file_to_counts(full_filename, bam=options.bam, extra=options.extra_out,
                                 use_all_references=options.use_all_references, min_mapq=options.min_mapq)

        if options.extra_out:
            counts, unique, nonunique = tmp['counts'], tmp['unique'], tmp['nonunique']
        else:
            counts = tmp['counts']

        ## save counts, and unique/non-unique counts
        file_counts[filename] = counts

        if options.extra_out:
            file_unique_counts[filename] = unique
            file_nonunique_counts[filename] = nonunique

        ## add all ids encountered in this in this file
        all_ids.extend(file_counts[filename].keys())

    ## Uniquify all_ids, and then take the nested file_counts
    ## dictionary, collapse, and write to file.
    all_ids = set(all_ids)
    table_dict = collapsed_nested_count_dict(file_counts, all_ids, order=files)
    counts_to_file(table_dict, options.out_file, delimiter=options.delimiter)

    if options.extra_out:
        unique_fn, nonunique_fn = options.extra_out_files.split(',')
        unique_table_dict = collapsed_nested_count_dict(file_unique_counts, all_ids, order=files)
        nonunique_table_dict = collapsed_nested_count_dict(file_nonunique_counts, all_ids, order=files)
        
        counts_to_file(unique_table_dict, unique_fn, delimiter=options.delimiter)
        counts_to_file(nonunique_table_dict, nonunique_fn, delimiter=options.delimiter)
        
