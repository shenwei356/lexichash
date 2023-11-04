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

package lexichash

import (
	"strconv"
	"testing"
	"time"

	"github.com/shenwei356/bio/seqio/fastx"
)

func TestIndex(t *testing.T) {
	k := 21
	nMasks := 1000
	idx, err := NewIndexWithSeed(k, nMasks, 11)
	if err != nil {
		t.Error(err)
		return
	}

	queries, err := fastx.GetSeqs("tests/hairpin.query.fasta", nil, 8, 100, "")
	if err != nil {
		t.Error(err)
		return
	}
	if len(queries) == 0 {
		t.Error(err)
		return
	}
	// queryID := queries[0].ID

	seqs, err := fastx.GetSeqs("tests/hairpin.fasta", nil, 8, 100, "")
	if err != nil {
		t.Error(err)
		return
	}

	sTime := time.Now()
	t.Logf("starting to build the index ...")
	for _, s := range seqs {
		// if bytes.Equal(s.ID, queryID) { //skip the query sequence
		// 	continue
		// }
		idx.Insert(s.ID, s.Seq.Seq)
	}
	t.Logf("finished to build the index in %s from %d sequences with %d masks",
		time.Since(sTime), len(seqs), nMasks)

	for _, s := range queries {
		sr, err := idx.Search(s.Seq.Seq, 12)
		if err != nil {
			t.Log(err)
			return
		}
		t.Log()
		t.Logf("query: %s, targets: %d\n", s.ID, len(sr))
		if sr == nil {
			continue
		}

		for i, r := range sr {
			t.Logf("%4s %s\n", "#"+strconv.Itoa(i+1), idx.IDs[r.IdIdx])
			for _, v := range r.Subs[0:] {
				t.Logf("     k:%d %s (%d-%d,%v) vs (%d-%d,%v)\n",
					v[0].K, v[0].KmerCode,
					v[0].Begin+1, v[0].End, v[0].RC,
					v[1].Begin+1, v[1].End, v[1].RC)
			}
		}
	}

}
