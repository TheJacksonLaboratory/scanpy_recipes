#!/usr/bin/env python
import os
import argparse
import nbformat

parser = argparse.ArgumentParser()
parser.add_argument("notebook_in")
parser.add_argument("notebook_out")
parser.add_argument("remove_tags", nargs="+")
args = parser.parse_args()

nbin = nbformat.read(args.notebook_in, nbformat.NO_CONVERT)
for cell in nbin.cells.copy():
    tags = cell.metadata.get("tags", [])
    for remove_tag in args.remove_tags:
        if remove_tag in tags:
            nbin.cells.remove(cell)
            break
nbformat.write(nbin, args.notebook_out)
