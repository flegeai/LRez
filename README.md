# LRez - playing with 10X fastq and bam files 

***A new and more efficient version of the tools is available at <http://github.com/morispi/LREz>***

[![License](http://img.shields.io/:license-affero-blue.svg)](http://www.gnu.org/licenses/agpl-3.0.en.html)

LRez proposes few tools for working with 10x barcoded reads :

1. *idx_bx_sqlite3.py* : indexation iby barcodes of fastq.gz files x(idx_bx_sqlite3.py)
2. *reads_bx_sqlite3.py* extraction of reads corresponding to a list of barcodes from indexed fastq files (with idx_bx_sqlite3.py)
3. *BamExtractor* : extracts barcode from regions of bam files 
4. *Compare* : counts the common barcodes of regions from bam files 


### Requirements

1. C++
htslib 
cmake (3.4+)

2. Python3
shelve
sqlite3
indexed_gzip

### Usage

#### idx_bx_sqlite3.py 
    python3 idx_bx_sqlite3.py [-h] -bx BASIC -idx IDX [-z] [-m MODE]
	BASIC : barcoded Fastq file from reads obtained with longranger basic 
	IDX : output indexed file
	MODE: shelve/sqlite mode of indexation (default : sqlite)
	-z (--gz) : the fastq is zipped (defauld true)

#### reads_bx_sqlite3.py
    python3 reads_bx_sqlite3.py [-h] -f FASTQ -i IDX -b BDX [-z] [-m MODE]
	FASTQ : indexed fastq file 
	IDX : index file generated by idx_bx_sqlite3.py
	BX : list of barcodes 
	MODE: shelve/sqlite mode of indexation (default : sqlite)
        -z (--gz) : the fastq is zipped (defauld true)

#### BamExtractor
    BamExtractor BAM [REGION]
	BAM : a samtools indexed  bam file
	REGION : a specific region (e.g. Scaffold123:100000-200000) 

#### Compare
    Compare --bam BAM [--list LIST] [--in : CONTIG] [--size BOUND]
	BAM : a samtools indexed  bam file
	LIST : a list of regions to be compared in a pairwise mode
	CONTIG : a contig name (to be compare to all the others targets of the bam file 
	BOUND : boundaries size of the CONTIG option (only boundaries region are taken einto account) default = 1000

### Installation

    git clone https://github.com/flegeai/LRez.git

LRez is also distributed as a [Bioconda package](https://anaconda.org/bioconda/lrez):

	conda install -c bioconda lrez	

### Contact

LRez is a tool developed by Fabrice Legeai <fabrice.legeai@inrae.fr>
