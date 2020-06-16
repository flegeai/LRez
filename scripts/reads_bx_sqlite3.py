#!/usr/bin/env python3

"""*******************************************************************************
    Name: LRez/reads_bx_sqlite3.py 
    Description: extracts reads with barcodes from a set of indexed fastq file with idx_bx_sqlite3.py 
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
import indexed_gzip as igzip
import re
import sqlite3
import shelve
import sys

# Construct the argument parser
ap = argparse.ArgumentParser()
# Add the arguments to the parser
ap.add_argument("-f", "--fastq", required=True, help="gzipped barcoded Fastq file from reads obtained with longranger basic")
ap.add_argument("-i", "--idx", required=True, help="barcodes index file")
ap.add_argument("-b", "--bdx",  required=True, help="barcodes list")
ap.add_argument("-z", "--gz", help="fastq is zipped", action='store_true')
ap.add_argument("-m", "--mode", help="mode of storage (shelve/sqlite)", default="sqlite" )
args = vars(ap.parse_args())

if (args["mode"]=="sqlite"):
    conn = sqlite3.connect(args["idx"])
    c = conn.cursor()
if (args["mode"]=="shelve"):
    pos_bar=shelve.open(args["idx"],writeback=False)

#myfile = open(args["fastq"], "rt")
if (args["gz"]):
    myfile = igzip.IndexedGzipFile(args["fastq"], index_file= args["fastq"]+'.gzidx')
else:
    myfile=open(args["fastq"], "rt")

bx_re=re.compile("BX:Z:([ATCG]{16})-1")
for barcode in open(args['bdx']):
    print("barcode "+ barcode, file=sys.stderr)
    barcode=barcode.rstrip()
    pos =-1
    if (args["mode"]=="sqlite"):
        t = (barcode,)
        c.execute('select pos from bx_pos where barcode=?', t)
        for row in c:
            pos=row[0]

    if (args["mode"]=="shelve"):
        if barcode in pos_bar:
            pos=pos_bar[barcode]

    compt=0
    got_them_all=0

    if (pos != -1):
        myfile.seek(pos)

        if (args["gz"]):
            line = myfile.readline().decode('UTF-8').rstrip()
        else:
            line = myfile.readline().rstrip()

        while line:
            if (compt%4==0):
                match=bx_re.search(line)
                if (match != None):
                    bx=match.group(1)
                    if (bx != barcode):
                        got_them_all=1
                        break
                else:
                    got_them_all=1
                    break

            print(line)

            if (args["gz"]):
                line = myfile.readline().decode('UTF-8').rstrip()
            else:
                line = myfile.readline().rstrip()
        #        print(line,file=sys.stderr)
            compt+=1

        if (got_them_all != 1):
            print("Missing reads from barcode " + barcode, file=sys.stderr)

if (args["mode"]=="sqlite"):
    c.close()
