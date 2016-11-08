#!/bin/bash
set -e;

echo "# CHAR-LEVEL SEGMENTATION INDEX:";
../kaldi-lattice-word-index --symbols-table=syms.txt 3 ark:lattice.txt

echo "# CHAR-LEVEL SEGMENTATION INDEX WITH BEST SEGMENTATION/WORD:";
../kaldi-lattice-word-index \
  --only-best-segmentation --symbols-table=syms.txt 3 ark:lattice.txt


echo "# WORD-LEVEL SEGMENTATION INDEX:"
../kaldi-lattice-word-index --symbols-table=syms.txt 3 ark:lattice2.txt |
awk '{
  # lattice id, log-prob, word length (in characters)
  lat=$1; logp=$2; len = (NF - 3) / 2;
  # start and end frames of the whole word
  f0 = $(3 + len); f1 = $NF;
  # Get word
  word = $3;
  for (i = 1; i < len; ++i) word=word" "$(i + 3);
  token = lat"____"word"____"f0"____"f1;
  if (token in IDX) {
    if (logp > IDX[token]) { tmp = logp; logp = IDX[token]; IDX[token] = tmp; }
    IDX[token] += log(1.0 + exp(logp - IDX[token]));
  } else {
    IDX[token] = logp;
  }
}END{
  for (token in IDX) {
    logp = IDX[token];
    split(token, arr, "____");
    print arr[1], logp, arr[2], arr[3], arr[4];
  }
}'
