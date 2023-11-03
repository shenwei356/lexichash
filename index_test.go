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
	"testing"
	"time"

	"github.com/shenwei356/bio/seqio/fastx"
)

func TestIndex(t *testing.T) {
	k := 31
	nMasks := 200
	idx, err := NewIndex(k, nMasks)
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
	// fmt.Printf("number of elements in trees:\n")
	// for i, tree := range idx.Trees {
	// 	fmt.Printf("tree #%d: %d\n", i, tree.Len())
	// }
	t.Logf("finished to build the index in %s", time.Since(sTime))

	for _, s := range queries {
		idx.Search(s.Seq.Seq, 10)
	}

}
