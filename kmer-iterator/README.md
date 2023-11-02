# kmer-iterator

This package is very similar to the [k-mer iterator](https://github.com/shenwei356/bio/blob/master/sketches/iterator.go).
But it only supports k<=31, with the last two bits as a flag indicating if the k-mer is from the negative strand.
