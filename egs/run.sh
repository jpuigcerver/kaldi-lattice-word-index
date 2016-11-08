#!/bin/bash
set -e;

echo "INDEX WITH MULTIPLE WORD OCCURRENCES/LATTICE:"
../kaldi-lattice-word-index --symbols-table=syms.txt 3 ark:lattice.txt

echo "INDEX WITH BEST-WORD OCCURRENCE/LATTICE:";
../kaldi-lattice-word-index \
  --only-best-segmentation --symbols-table=syms.txt 3 ark:lattice.txt
