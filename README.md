# kaldi-lattice-word-index

This tool builds a word index from character lattices.

Build a word index from character lattices.

Words are any sequence of characters in between any of the separator symbols
(e.g.: whitespace, punctuation marks, etc).

The program will output the n-best character segmentations of words, with their
scores. More precisely:

Let's define a binary variable R_c that denotes whether the transcription (y) of
a sample (x) contains the word formed by characters c_{1:n}.

![R_c = 1 \iff y \in (\Sigma^* S)? c_{1:n} (S \Sigma^*)?](http://www.sciweavers.org/tex2img.php?eq=R_c%20%3D%20%5Cbegin%7Bcases%7D%0A1%20%26%20%5Cmathbf%7By%7D%20%5Cin%20%28%5CSigma%5E%2A%20%7E%20S%29%5E%3F%20%5Cmathbf%7Bc%7D%20%28S%20%7E%20%5CSigma%5E%2A%29%5E%3F%5C%5C%0A0%20%26%20%5Ctext%7Botherwise%7D%0A%5Cend%7Bcases%7D&bc=White&fc=Black&im=jpg&fs=12&ff=modern&edit=0)

Let s_{1:n} be a segmentation of each character in c_{1:n}, then the program
computes:

- If --only-best-segmentation=false (the default) then:

![n-best_{c_{1:n},s_{1:n}} P(R_c = 1| x, s_{1:n})](http://www.sciweavers.org/tex2img.php?eq=%5Ctext%7Bn-best%7D_%7B%5Cmathbf%7Bc%7D%2C%20%5Cmathbf%7Bs%7D%7D%20P%28R_c%20%3D%201%20%5Cmid%20%5Cmathbf%7Bx%7D%2C%20%5Cmathbf%7Bs%7D%29&bc=White&fc=Black&im=jpg&fs=12&ff=modern&edit=0)

- If --only-best-segmentation=true then:

![n-best_{c_{1:n},s_{1:n}} max_{s_{1:n}} P(R_c = 1 | x, s_{1:n})](http://www.sciweavers.org/tex2img.php?eq=%5Ctext%7Bn-best%7D_%7B%5Cmathbf%7Bc%7D%2C%20%5Cmathbf%7Bs%7D%7D%20%5Cmax_%7B%5Cmathbf%7Bs%7D%7D%20P%28R_c%20%3D%201%20%5Cmid%20%5Cmathbf%7Bx%7D%2C%20%5Cmathbf%7Bs%7D%29&bc=White&fc=Black&im=jpg&fs=12&ff=modern&edit=0)

This gives a lower bound to P(R_c = 1 | x), but it is usually quite close.


## Usage
```
Usage: kaldi-lattice-word-index [options] separator-symbols lat-rspecifier
 e.g.: kaldi-lattice-word-index "1 2" ark:lats.ark
 e.g.: kaldi-lattice-word-index --nbest=10000 "1 2" ark:lats.ark
`````

## Options
  - --acoustic-scale            : Scaling factor for acoustic likelihoods in the lattices. (float, default = 1)
  - --beam                      : Pruning beam (applied after acoustic scaling and adding the insertion penalty). (float, default = inf)
  - --delta                     : Tolerance used in determinization. (float, default = 0.000976562)
  - --graph-scale               : Scaling factor for graph probabilities in the lattices. (float, default = 1)
  - --insertion-penalty         : Add this penalty to the lattice arcs with non-epsilon output label (typically, equivalent to word insertion penalty). (float, default = 0)
  - --max-mem                   : Maximum approximate memory usage in determinization (real usage might be many times this). (int, default = 536870912)
  - --nbest                     : Extract this number of n-best hypothesis. (int, default = 100)
  - --only-best-segmentation    : If true, output the best character segmentation for each word. (bool, default = false)
  - --symbols-table             : Use this symbols table to map from labels to characters. (string, default = "")
