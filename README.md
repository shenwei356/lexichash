# lexichash

[![Go Reference](https://pkg.go.dev/badge/github.com/shenwei356/lexichash.svg)](https://pkg.go.dev/github.com/shenwei356/lexichash)

This project implements [LexicHash](https://academic.oup.com/bioinformatics/article/39/11/btad652/7329717) in Golang.

## Command-line tool

A simple tool `lexic-hash-search` is available for finding shared substrings between query and target sequences.

Example:

```
$ memusg -t -s "lexic-hash-search -k 21 -n 1000 -m 13 -j 16 tests/hairpin.query.fasta tests/hairpin.fasta \
    | csvtk pretty -t"
2023/11/04 23:38:53 starting to build the index from 1 files
2023/11/04 23:39:08 finished building the index in 15.031332246s from 26379 sequences with 1000 masks
2023/11/04 23:39:08 finished searching with 1 sequences in 2.254832ms
query   target         qstart   qend   qstrand   tstart   tend   tstrand   len   match
-----   ------------   ------   ----   -------   ------   ----   -------   ---   --------------------
query   cre-MIR9897    1        20     +         1        20     +         20    AGAAGGACGTGGACGTGGAT
query   cre-MIR9897    22       34     +         22       34     +         13    CCGATAAGAAGGA
query   cre-MIR9897    36       48     -         35       47     -         13    GGTACCTTACGGC
query   cre-MIR9897    50       68     -         49       67     -         19    CCCCTGCCCTCCCCACGCC
query   cre-MIR9897    70       88     -         69       87     -         19    GCCCCTGATCCCCGTCCCT
query   cre-MIR9897    71       90     +         70       89     +         20    GGGACGGGGATCAGGGGCAG
query   ptr-mir-3138   57       72     -         7        22     +         16    CCTTCCCCTGCCCTCC
query   gma-MIR5376    28       40     -         22       34     -         13    ACGGCGTCCTTCT
query   hsa-mir-6850   64       76     +         25       37     +         13    AGGGGAAGGGACG
query   hsa-mir-6799   59       71     -         40       52     +         13    CTTCCCCTGCCCT
query   bta-mir-7865   57       69     +         4        16     +         13    GGAGGGCAGGGGA
query   gga-mir-1456   65       77     +         74       86     -         13    GGGGAAGGGACGG


elapsed time: 15.107s
peak rss: 3.32 GB

$ seqkit stats tests/hairpin.fasta
file                 format  type  num_seqs    sum_len  min_len  avg_len  max_len
tests/hairpin.fasta  FASTA   DNA     26,379  2,748,673       39    104.2    2,354
```

Pairwise searching

```
$ memusg -t -s "lexic-hash-search -k 21 -n 1000 -m 15 -j 16 tests/hairpin.fasta tests/hairpin.fasta | pigz -c > t.gz"
2023/11/04 23:45:22 starting to build the index from 1 files
2023/11/04 23:45:36 finished building the index in 13.583252685s from 26379 sequences with 1000 masks
2023/11/04 23:46:39 finished searching with 26379 sequences in 1m3.360959079s

elapsed time: 1m:17s
peak rss: 7.65 GB
```

Usage

```
This command uses LexicHash to find shared substrings between query and target sequences.

Author: Wei Shen <shenwei356@gmail.com>
  Code: https://github.com/shenwei356/lexichash

Version: v0.1.0
Usage: lexic-hash-search [options] <query fasta/q> <target fasta/q> [<target fasta/q> ...]

Options/Flags:
  -c    using cannocial k-mers
  -h    print help message
  -j int
        number of threads (default 16)
  -k int
        k-mer size (default 21)
  -m int
        minimum length of shared substrings (default 13)
  -n int
        number of maskes/hashes (default 1000)
  -s int
        seed number (default 1)
```


## Support

Please [open an issue](https://github.com/shenwei356/lexichash/issues) to report bugs,
propose new functions or ask for help.

## License

[MIT License](https://github.com/shenwei356/lexichash/blob/master/LICENSE)

