// Copyright © 2023-2024 Wei Shen <shenwei356@gmail.com>
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//b
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.

package main

import (
	"errors"
	"flag"
	"fmt"
	"io"
	"log"
	"os"
	"path/filepath"
	"runtime"
	"time"

	"github.com/pkg/profile"
	"github.com/shenwei356/bio/seq"
	"github.com/shenwei356/bio/seqio/fastx"
	"github.com/shenwei356/lexichash"
	"github.com/shenwei356/xopen"
)

var version = "0.1.0"

func main() {
	usage := fmt.Sprintf(`
This command uses LexicHash to find shared substrings between query and target sequences.

Author: Wei Shen <shenwei356@gmail.com>
  Code: https://github.com/shenwei356/lexichash

Version: v%s
Usage: %s [options] <query fasta/q> <target fasta/q> [<target fasta/q> ...]

Options/Flags:
`, version, filepath.Base(os.Args[0]))

	flag.Usage = func() {
		fmt.Fprint(os.Stderr, usage)
		flag.PrintDefaults()
	}

	help := flag.Bool("h", false, "print help message")
	k := flag.Int("k", 21, "k-mer size")
	nMasks := flag.Int("n", 1000, "number of maskes/hashes")
	seed := flag.Int("s", 1, "seed number")
	cannonical := flag.Bool("c", false, "using cannocial k-mers")
	minLen := flag.Int("m", 13, "minimum length of shared substrings")
	threads := flag.Int("j", runtime.NumCPU(), "number of threads")
	pfCPU := flag.Bool("pprof-cpu", false, "pprofile CPU")
	pfMEM := flag.Bool("pprof-mem", false, "pprofile memory")

	flag.Parse()

	if *help {
		flag.Usage()
		return
	}

	if flag.NArg() < 2 {
		flag.Usage()
		os.Exit(1)
	}

	if *k > 31 {
		checkError(fmt.Errorf("k should be < 31"))
	}
	if *nMasks < 1 {
		checkError(fmt.Errorf("n should be >=1, >=100 or >=1000 recommended"))
	}
	if *seed < 1 {
		checkError(fmt.Errorf("s should be > 0"))
	}
	if *minLen > *k {
		checkError(fmt.Errorf("m should be <= k"))
	}
	if *threads <= 0 {
		*threads = runtime.NumCPU()
	}

	for _, file := range flag.Args() {
		if _, err := os.Stat(file); errors.Is(err, os.ErrNotExist) {
			checkError(fmt.Errorf("%s", err))
		}
	}

	// -----------------------------------------------

	// go tool pprof -http=:8080 cpu.pprof
	if *pfCPU {
		defer profile.Start(profile.CPUProfile, profile.ProfilePath(".")).Stop()
	} else if *pfMEM {
		defer profile.Start(profile.MemProfile, profile.ProfilePath(".")).Stop()
	}

	outfh, err := xopen.Wopen("-")
	checkError(err)
	defer outfh.Close()

	lexichash.Threads = *threads
	idx, err := lexichash.NewIndexWithSeed(*k, *nMasks, *cannonical, int64(*seed))
	checkError(err)

	log.Printf("starting to build the index from %d files", len(flag.Args()[1:]))
	sTime := time.Now()

	seq.ValidateSeq = false
	var record *fastx.Record
	var fastxReader *fastx.Reader
	var nSeqs int
	for _, file := range flag.Args()[1:] {
		fastxReader, err = fastx.NewReader(nil, file, "")
		checkError(err)

		for {
			record, err = fastxReader.Read()
			if err != nil {
				if err == io.EOF {
					break
				}
				checkError(err)
				break
			}

			nSeqs++

			idx.Insert([]byte(string(record.ID)), record.Seq.Seq)
		}
	}
	log.Printf("finished building the index in %s from %d sequences with %d masks",
		time.Since(sTime), nSeqs, *nMasks)

	// -----------------------------------------------

	sTime = time.Now()

	fmt.Fprintf(outfh, "query\ttarget\tqstart\tqend\tqstrand\ttstart\ttend\ttstrand\tlen\tmatch\n")

	fastxReader, err = fastx.NewReader(nil, flag.Args()[0], "")
	checkError(err)
	var sr *[]*lexichash.SearchResult
	var r *lexichash.SearchResult
	var v *[2]lexichash.Substr
	var nQueries int
	decoder := lexichash.MustDecoder()

	for {
		record, err = fastxReader.Read()
		if err != nil {
			if err == io.EOF {
				break
			}
			checkError(err)
			break
		}

		nQueries++

		sr, err = idx.Search(record.Seq.Seq, uint8(*minLen))
		if err != nil {
			checkError(err)
		}
		if sr == nil {
			continue
		}

		for _, r = range *sr {
			for _, v = range *r.Subs {
				fmt.Fprintf(outfh, "%s\t%s\t%d\t%d\t%c\t%d\t%d\t%c\t%d\t%s\n",
					record.ID, idx.IDs[r.IdIdx],
					v[0].Begin+1, v[0].End, lexichash.Strands[v[0].RC],
					v[1].Begin+1, v[1].End, lexichash.Strands[v[1].RC],
					v[0].K, decoder(v[0].Code, v[0].K))
			}
		}
		idx.RecycleSearchResult(sr)
	}

	log.Printf("finished searching with %d sequences in %s",
		nQueries, time.Since(sTime))
}

func checkError(err error) {
	if err != nil {
		fmt.Println(err)
		os.Exit(1)
	}
}
