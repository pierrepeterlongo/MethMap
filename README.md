

MethMap
========
[![License](http://img.shields.io/:license-affero-blue.svg)](http://www.gnu.org/licenses/agpl-3.0.en.html)
[![Build Status](https://travis-ci.org/pierrepeterlongo/MethMap.svg?branch=master)](https://travis-ci.org/pierrepeterlongo/MethMap)
# What is MethMap?

MethMap is designed for mapping short sequences (ie reads) to a reference bank. If required, when mapping converted sequences, 'T's from queries match 'C's from the bank. Otherwise, this script is a classical seed-and-extend heuristic.

# Getting the latest source code
**Requirements**

> python 3

**Instructions**

get a local copy of MethMap source code
> git clone --recursive https://github.com/pierrepeterlongo/MethMap.git

Run a simple test on your computer
> cd MethMap
> cd tests
> ./test.sh
> cd ..


# User manual

## MethMap usage:
> MethMap.py [-h] [-k K] [-t T] [-span S] input_bank_file input_query_file converted

- positional arguments:
    - input_bank_file   input fasta or fastq bank file
    - input_query_file  input fasta or fastq query file
    - converted         convert: chose "True" or "False". False: usual mapping.
                        True: "T"s from queries match "C"s from the bank.
- optional arguments:
    - -h, --help        show this help message and exit
    - -k K              kmer size [default: 12]
    - -t T              Maximal number authorized substitution [default: 0]
    - -span S           The portion of a read mapped on a reference may be lower
                        than 100 percent. Span (in 0-100) provides this minimal
                        percentage value [Default 90]




####Warning
 No reverse complement is considered. All input sequences are considered 5' 3'

####Contacts
 - Amandine Etcheverry:  amandine.etcheverry@univ-rennes1.fr
 - Marc Aubry: marc.aubry@univ-rennes1.fr
 - Pierre Peterlongo: pierre.peterlongo@inria.fr