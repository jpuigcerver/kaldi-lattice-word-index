# kaldi-lattice-word-index

This tool builds a word index from character lattices.

Words are any sequence of symbols in between two separators (typical cases: whitespace, punctuation marks, etc).

The scores of each word are lower bounds of the probability that the word 
is present in the transcription (according to the lattice hypotheses), but 
typically they are very close to the exact probability.

## Usage
```
Usage: kaldi-lattice-word-index [options] separator-symbols lat-rspecifier
 e.g.: kaldi-lattice-word-index "1 2" ark:lats.ark
 e.g.: kaldi-lattice-word-index --nbest=10000 "1 2" ark:lats.ark
`````

## Options
  - --acoustic-scale       : Scaling factor for acoustic likelihoods in the lattices. (float, default = 1)
  - --beam                 : Pruning beam (applied after acoustic scaling and adding the insertion penalty). (float, default = inf)
  - --determinize          : Determinize fst before building the index. (bool, default = true)
  - --graph-scale          : Scaling factor for graph probabilities in the lattices. (float, default = 1)
  - --insertion-penalty    : Add this penalty to the lattice arcs with non-epsilon output label (typically, equivalent to word insertion penalty). (float, default = 0)
  - --nbest                : Extract this number of paths from the words fst. If the word fst is deterministic (--determinize=true), then this is the number of words in the index. Otherwise, the actual number of words in the index can be smaller. (int, default = 100)
