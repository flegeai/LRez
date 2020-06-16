#!/usr/bin/env python3

"""*******************************************************************************
    Name: LRez/idx_bx_sqlite3.py
    Description: index 10x reads sets (in fastq.gz) by barcodes   
    Version: 1.0.0
    Author: Fabrice Legeai
    Contact: fabrice.legeai@inrae.fr, IGEPP, INRAE, Institut Agro, Univ Rennes, 35000, Rennes, France
    Copyright (C) 2020 INRAE
    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU Affero General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Affero General Public License for more details.
    You should have received a copy of the GNU Affero General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
*******************************************************************************"""
__version__ = '1.0.0'

import argparse
import gzip
import re
import sqlite3
import shelve
import indexed_gzip as igzip

# Construct the argument parser
ap = argparse.ArgumentParser()
# Add the arguments to the parser
ap.add_argument("-bx", "--basic", required=True, help="barcoded Fastq file from reads obtained with longranger basic")
ap.add_argument("-idx", "--idx", required=True, help="output indexed file")
ap.add_argument("-z", "--gz", help="fastq is zipped", action='store_true')
ap.add_argument("-m", "--mode", help="mode of storage (shelve/sqlite)", default="sqlite" )
args = vars(ap.parse_args())

# First we are generating an index file for best performance


# Load the file, pre-generate the
# index, and save it out to disk.
if (args["gz"]):
    fobj = igzip.IndexedGzipFile(args["basic"])
    fobj.build_full_index()
    fobj.export_index(args["basic"]+'.gzidx')

# Get the barcode positions
if (args["mode"]=="sqlite"):
    conn = sqlite3.connect(args["idx"])
    c = conn.cursor()
    c.execute('create table bx_pos (barcode text, pos int)')
if (args["mode"]=="shelve"):
    pos_bar=shelve.open(args["idx"],writeback=True)

current_bx=''
compt=0
bx_re=re.compile("BX:Z:([ATCG]{16})-1")
if (args["gz"]):
    f=gzip.open(args["basic"], "rt")
else:
    f=open(args["basic"], "rt")

pos=f.tell()
line = f.readline()
while line:
    if (compt%4==0):
        match=bx_re.search(line)
        if (match != None):
            bx=match.group(1)
            if (bx != current_bx):
                #print(bx + " " + str(pos))
                # Insert a row of data
                if (args["mode"]=="sqlite"):
                    t=(bx, pos)
                    c.execute('insert into bx_pos values (?,?)', t)
                if (args["mode"]=="shelve"):
                    pos_bar[bx]=pos
                current_bx=bx
    pos=f.tell()
    line=f.readline()
    compt+=1

if (args["mode"]=="sqlite"):
    conn.commit()
    c.close()
if (args["mode"]=="shelve"):
    pos_bar.close()
