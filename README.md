sam2counts.py - Convert SAM mapping results to reference sequence counts
========================================================================

Requirements
------------
 - Python 2.6.5 or greater
 - Cython (download from <http://cython.org/>)
 - Pysam (download from <http://code.google.com/p/pysam/>)

Usage
-----

Bare minimum: 

    python sam2counts.py sample_1.sam sample_2.sam

This produces a tab-delimited file named counts.txt that contains the
counts for each reference. If you prefer to have it named something else:

    python sam2counts.py -o s1_to_s2.txt sample_1.sam sample_2.sam

Or use a different delimiter:

    python sam2counts.py -o s1_to_s2.csv -d, sample_1.sam sample_2.sam

Todo
-----

Currently there is code to indicate if mapped reads map with 0 quality
(i.e. map to multiple locations). The best output for this is still
being decided, and this is not entirey implemented

Furthermore, this code was tested fairly heavily, but I make no
guarantees.

License and Contact
----------------

If you have questions regarding this code, please contact Vince
Buffalo at vsbuffaloAAAAA@gmail.com (with the poly-A tail
removed). All code copyright Vince Buffalo, 2010.

License: GNU General Public License
