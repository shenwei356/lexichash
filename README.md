# lexichash

[![Go Reference](https://pkg.go.dev/badge/github.com/shenwei356/lexichash.svg)](https://pkg.go.dev/github.com/shenwei356/lexichash)

This project implements [LexicHash](https://academic.oup.com/bioinformatics/article/39/11/btad652/7329717) in Golang.

## Command-line tool

A simple tool `lexic-hash-search` is available for finding shared substrings between query and target sequences.

Example:

```
$ memusg -t -s "lexic-hash-search -k 21 -n 1000 -m 13 -j 16 tests/hairpin.query.fasta tests/hairpin.fasta \
    | csvtk pretty -t"
2023/11/06 20:40:23 starting to build the index from 1 files
2023/11/06 20:40:36 finished building the index in 13.259185495s from 26379 sequences with 1000 masks
2023/11/06 20:40:36 finished searching with 1 sequences in 224.866024ms
query   target         qstart   qend   qstrand   tstart   tend   tstrand   len   match
-----   ------------   ------   ----   -------   ------   ----   -------   ---   --------------------
query   cre-MIR9897    1        20     +         1        20     +         20    AGAAGGACGTGGACGTGGAT
query   cre-MIR9897    22       34     +         22       34     +         13    CCGATAAGAAGGA
query   cre-MIR9897    36       48     +         126      138    -         13    GCCGTAAGGTACC
query   cre-MIR9897    50       68     -         49       67     -         19    CCCCTGCCCTCCCCACGCC
query   cre-MIR9897    70       88     -         69       87     -         19    GCCCCTGATCCCCGTCCCT
query   cre-MIR9897    71       90     +         70       89     +         20    GGGACGGGGATCAGGGGCAG
query   ptr-mir-3138   57       72     +         7        22     -         16    GGAGGGCAGGGGAAGG
query   gga-mir-1456   65       77     +         74       86     -         13    GGGGAAGGGACGG
query   gma-MIR5376    28       40     -         22       34     -         13    ACGGCGTCCTTCT
query   hsa-mir-6799   59       71     -         40       52     +         13    CTTCCCCTGCCCT
query   hsa-mir-6850   64       76     +         25       37     +         13    AGGGGAAGGGACG
query   bta-mir-7865   57       69     +         4        16     +         13    GGAGGGCAGGGGA

elapsed time: 13.620s
peak rss: 3.36 GB

$ seqkit stats tests/hairpin.fasta
file                 format  type  num_seqs    sum_len  min_len  avg_len  max_len
tests/hairpin.fasta  FASTA   DNA     26,379  2,748,673       39    104.2    2,354
```

Pairwise searching

```
$ memusg -t -s "lexic-hash-search -k 21 -n 1000 -m 15 -j 16 tests/hairpin.fasta tests/hairpin.fasta | pigz -c > t.gz"
2023/11/06 20:41:41 starting to build the index from 1 files
2023/11/06 20:41:56 finished building the index in 14.969412816s from 26379 sequences with 1000 masks
2023/11/06 20:42:05 finished searching with 26379 sequences in 9.361483112s

elapsed time: 24.496s
peak rss: 4.13 GB
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
  -pprof-cpu
        pprofile CPU
  -pprof-mem
        pprofile memory
  -s int
        seed number (default 1)
```


## Support

Please [open an issue](https://github.com/shenwei356/lexichash/issues) to report bugs,
propose new functions or ask for help.

## License

[MIT License](https://github.com/shenwei356/lexichash/blob/master/LICENSE)

