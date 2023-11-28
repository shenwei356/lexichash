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
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.

package index

import (
	"io"
	"strconv"
	"strings"
	"testing"
	"time"
	"unsafe"

	"github.com/shenwei356/bio/seq"
	"github.com/shenwei356/bio/seqio/fastx"
	"github.com/shenwei356/kmers"
	"github.com/shenwei356/lexichash"
)

func TestStructSize(t *testing.T) {
	t.Logf("struct: Sizeof, Alignof\n")
	t.Logf("SubstrPair: %d, %d", unsafe.Sizeof(SubstrPair{}), unsafe.Alignof(SubstrPair{}))
	t.Logf("SearchResult: %d, %d", unsafe.Sizeof(SearchResult{}), unsafe.Alignof(SearchResult{}))
}

func TestHash(t *testing.T) {
	k := 21
	nMasks := 1000
	var seed int64 = 1

	s1 := []byte("AGAAGGACGTGGACGTGGATGCCGATAAGAAGGAGCCGTAAGGTACCGGGCGTGGGGAGGGCAGGGGCAGGGACGGGGATCAGGGGCAGCTGATCCCCGT")
	s2 := []byte("AGAAGGACGTGGACGTGGATcCCGATAAGAAGGAcGCCGTAAGGTACCaGGCGTGGGGAGGGCAGGGGaAGGGACGGGGATCAGGGGCAGaTGATCCCCGT")

	minLen := 13

	// use the same sequence to build the index

	idx, err := NewIndexWithSeed(k, nMasks, seed)
	if err != nil {
		t.Error(err)
		return
	}
	idx.Insert([]byte("s1"), s1, len(s1), []int{len(s1)})

	for i, info := range idx.RefSeqInfos {
		t.Logf("%s; sum: %d, concatenated sum: %d, #seqs: %d, #sizes: %v",
			idx.IDs[i], info.GenomeSize, info.Len, info.NumSeqs, info.SeqSizes)
	}

	sr, err := idx.Search(s2, uint8(minLen))
	if err != nil {
		t.Log(err)
		return
	}
	t.Log()
	t.Logf("query: %s", "s2")
	for i, r := range *sr {
		t.Logf("%4s %s\n", "#"+strconv.Itoa(i+1), idx.IDs[r.IdIdx])
		for _, v := range *r.Subs {
			t.Logf("     (%3d,%3d, %c) vs (%3d,%3d, %c) %3d %s\n",
				v.QBegin+1, v.QEnd, Strands[v.QRC],
				v.TBegin+1, v.TEnd, Strands[v.TRC],
				v.QK, kmers.MustDecode(v.QCode, int(v.QK)))
		}
	}
	idx.RecycleSearchResult(sr)
}

func TestIndex(t *testing.T) {
	k := 21
	nMasks := 1000
	var seed int64 = 1

	idx, err := NewIndexWithSeed(k, nMasks, seed)
	if err != nil {
		t.Error(err)
		return
	}

	queries, err := fastx.GetSeqs("../test_data/hairpin.query.fasta", nil, 8, 100, "")
	if err != nil {
		t.Error(err)
		return
	}
	if len(queries) == 0 {
		t.Error(err)
		return
	}

	sTime := time.Now()
	t.Logf("starting to build the index ...")

	input, done := idx.BatchInsert()

	seq.ValidateSeq = false
	var record *fastx.Record
	fastxReader, err := fastx.NewReader(nil, "../test_data/hairpin.fasta", "")
	if err != nil {
		t.Error(err)
		return
	}
	var nSeqs int
	for {
		record, err = fastxReader.Read()
		if err != nil {
			if err == io.EOF {
				break
			}
			t.Error(err)
			return
		}

		nSeqs++

		_seq := make([]byte, len(record.Seq.Seq))
		copy(_seq, record.Seq.Seq)
		input <- RefSeq{
			ID:  []byte(string(record.ID)),
			Seq: _seq,
		}
	}
	close(input) // wait BatchInsert
	<-done       // wait BatchInsert

	t.Logf("finished to build the index in %s from %d sequences with %d masks",
		time.Since(sTime), nSeqs, nMasks)

	minLen := 13

	decoder := lexichash.MustDecoder()
	for _, s := range queries {
		sr, err := idx.Search(s.Seq.Seq, uint8(minLen))
		if err != nil {
			t.Log(err)
			return
		}
		if sr == nil {
			continue
		}
		t.Log()
		t.Logf("query: %s, targets: %d\n", s.ID, len(*sr))

		for i, r := range *sr {
			t.Logf("%4s %s\n", "#"+strconv.Itoa(i+1), idx.IDs[r.IdIdx])
			for _, v := range *r.Subs {
				t.Logf("     (%3d,%3d, %c) vs (%3d,%3d, %c) %3d %s\n",
					v.QBegin+1, v.QEnd, Strands[v.QRC],
					v.TBegin+1, v.TEnd, Strands[v.TRC],
					v.QK, decoder(v.QCode, v.QK))
			}
		}
		idx.RecycleSearchResult(sr)
	}

	// idx.Trees[666].Walk(func(code uint64, v []uint64) bool {
	// 	fmt.Printf("%s\n", kmers.MustDecode(code, k))
	// 	return false
	// })

	_queries := []string{
		"ACGGCTGGGAGATGGAGCCAG",
		"GCACATATACTACACACACAT",
	}
	for _, query := range _queries {
		code, _ := kmers.Encode([]byte(query))
		t.Log()
		t.Logf("path of %s\n", query)
		paths := idx.Paths(code, uint8(len(query)), uint8(len(query)))
		for _, path := range *paths {
			t.Logf("  tree: %d, prefix: %d, path: %s\n", path.TreeIdx, path.Bases, strings.Join(*path.Nodes, "->"))
		}
		idx.RecyclePathResult(paths)
	}
}
