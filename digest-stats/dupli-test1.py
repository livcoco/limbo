#!/usr/bin/env python

# Re: Duplication-testing for several file-extract sizes, using
# a block-digest method from hashlib module. - jiw - April 2018

# Note: To convert a program's python-2-style print statements to
# python-3-style, in place, you can use a command of following form:
#         2to3 -n -w <programname>

# This program reads and digests blocks of data from fronts of files
# under a given directory.  It tabulates and reports the numbers of
# duplications that occur when the digest numbers are truncated to
# different lengths.  For example, using python's sha1() hash
# algorithm (which produces 160-bit or L = 20-byte digests) as a
# digester, we report the number of collisions at lengths k = 1, 2,
# ... m bytes, [where m is either L or the least k with no collisions]
# and also report min and max bucket occupancy counts.  (At each
# length k, files whose first k digest bytes match are counted into
# the same bucket.)

# Optional command-line parameters:

#   BlockSize -- Max number of bytes to read from front of each
#       file; default 1024
#   fileDir -- Path to directory from which files are obtained.
#       Default is . (current working directory).
#       Note, we add a / at end of path if it's missing.
#       Note, files from subdirectories (as well as top level) will
#          be processed, if they exist and maxFiles is large.
#   maxFiles -- Max number of files to process.
#       Defaults to a small number.
#   snDig -- Show snDig first and last digest values.  Default 0

import time
from os  import walk
from sys import argv, stdout
import hashlib

# Program-specified parameter:
#   digester -- Method to use as digester
digester = hashlib.sha1
dibytes = digester().digest_size
diname = digester().name

# Get params:
arn = 0
arn+=1; blkSiz = int(argv[arn]) if len(argv)>arn else 1024
arn+=1; filDir = argv[arn]      if len(argv)>arn else '.'
arn+=1; maxFil = int(argv[arn]) if len(argv)>arn else 100
arn+=1; snDig  = int(argv[arn]) if len(argv)>arn else 0
if filDir[:-1] != '/':        # Need / at end for walk() to work right
    filDir = filDir + '/'

print 'Parameter values: blokSiz={}, fileDir="{}", maxFiles={}'.format(blkSiz, filDir, maxFil)
t0 = time.time()

# Make files list
w = walk(filDir)
flist = []
for dirpath, subdirnames, filenames in w:
    fpart = filenames[:maxFil-len(flist)]
    for f in fpart:
        flist.append(dirpath+'/'+f)
    if len(flist) >= maxFil:
        break

te = time.time()-t0
nfiles = len(flist)
print 'Listed   {} files in {:8.6f} seconds'.format(nfiles, te)
digs = []                       # Empty list of digests-so-far
ts = time.time()                # Record starting time

for fn in flist:
    with open(fn, 'r') as f:
        datablok = f.read(blkSiz)
        # Append digest bytes to digs array
        digs.append(digester(datablok).digest())
tr = time.time()-ts

print 'Computed {}  {}-byte {} digests in {:7.6f} seconds, or {:7.9f} sec. each'.format(nfiles, dibytes, diname, tr, tr/nfiles)

print
digs = sorted(digs)
if snDig:       # Optional print of some digest values
    print 'First {} and last {} digests:'.format(snDig, snDig)
    for i in range(snDig+1) + range(-snDig,0,1):
        d = digs[i]
        print ''.join('{:02x}'.format(ord(x)) for x in d)
    print

fmt = '{:>4} {:>6} {:>6} {:>6}'
print 'Comparing first B bytes of digests, putting B-byte matches'
print 'in buckets; with total matches and min/max bucket counts:'
print fmt.format('B', '#Dups', 'bMin', 'bMax')
for b in range(1, dibytes+1):
    dups = 0
    bmin = nfiles
    bmax = 0
    prev = digs[-1]
    inBuc = None
    for i in range(nfiles):
        if prev[0:b] == digs[i][0:b]:
            inBuc += 1
            dups  += 1
        else:
            if inBuc:
                bmax = max(bmax, inBuc)
                bmin = min(bmin, inBuc)
            inBuc = 1
        prev = digs[i]
    print fmt.format(b, dups, bmin, bmax)
    if dups==0:  break
print
