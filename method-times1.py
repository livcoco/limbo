#!/usr/bin/env python
# Re: Timing-tests of several block-digest methods from hashlib module
# - jiw - April 2018

# Note: To convert this program's python-2-style print
# statements to python-3-style, you can use the command
#         2to3 -w method-times1.py

# Optional parameters:

#   BlockSize -- Max number of bytes to read from each file; default 1024

#   passes  -- Number of times to time processing of the given set of files.
#              Default: 3.  (Extra passes may avoid or reveal initial costs)

#   fileDir -- Path to directory from which files are obtained.
#       Note, we add a / at end of path if it's missing.
#       Note, files from subdirectories (as well as top level) may
#          be processed, if they exist and maxFiles is large.

#   maxFiles -- Max number of files to process.  Defaults to a small number.

from sys import argv, stdout
import time
from hashlib import md5, sha1, sha224, sha256, sha384, sha512
from os  import walk
# Get params:
arn = 0
arn+=1; blkSiz = int(argv[arn]) if len(argv)>arn else 1024
arn+=1; passes = int(argv[arn]) if len(argv)>arn else 3
arn+=1; filDir = argv[arn]      if len(argv)>arn else '.'
arn+=1; maxFil = int(argv[arn]) if len(argv)>arn else 100
if filDir[:-1] != '/':
    filDir = filDir + '/'

print 'Parameter values: blokSiz={}, passes={}, fileDir="{}", maxFiles={}'.format(blkSiz, passes, filDir, maxFil)
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
print 'len flist = {} at {:8.6f} elapsed seconds'.format(nfiles, te)


def runtest(method):
    d = method()
    ts = time.time()            # Record starting time
    for fn in flist:
        with open(fn, 'r') as f:
            datablok = f.read(blkSiz)
            d = method(datablok)
    tr = time.time()-ts
    return tr, d.digest_size, d.name


print 'Pass -1 to pre-access files',
stdout.flush()
rtime, dbits, dname = runtest(md5)
print 'took {:8.6f} seconds'.format(rtime)
methodList = [md5, sha1, sha224, sha256, sha384, sha512]
ttimes = [0]*len(methodList)

for passnum in range(passes):
    print 'Pass {}:'.format(passnum),
    stdout.flush()
    if passnum<1:
        print
    for methNum, meth in enumerate(methodList):
        rtime, dbits, dname = runtest(meth)
        if passnum<1:
            print '{:2}.  {:3}-bit {:9} test: {:11.6f} seconds, {:11.9f} per file'.format(methNum, dbits, dname, rtime, rtime/nfiles)
        else:
            print dname,
            stdout.flush()
        ttimes[methNum] += rtime
    print
te = time.time()-t0
tfiles = nfiles*passes
print 'At end: {:8.6f} seconds total'.format(te)
print '\nTotal times and averages per file for each method, on {} files:'.format(tfiles)
for methNum, meth in enumerate(methodList):
    et = ttimes[methNum]
    d = meth()
    print '{:2}.  {:3}-bit {:9} test: {:11.6f} seconds, {:11.9f} per file'.format(methNum, d.digest_size, d.name, et, et/tfiles)
print; print
