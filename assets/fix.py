#!/usr/bin/env python
from os import rename, listdir
from sys import argv

badprefix = "cheese_"
fnames = listdir(argv[1])

for fname in fnames:
    if ".jpg" in fname:
        rename(fname,fname.replace('image_', ''))

print "Done"