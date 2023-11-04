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

package tree

import (
	"strings"
	"testing"

	"github.com/shenwei356/kmers"
)

func TestTree(t *testing.T) {
	tree := New()

	k := 4
	n := uint64(1 << (k * 2))
	var i uint64

	var v uint64
	var r []uint64
	var ok bool

	for i = 0; i < n; i++ {
		v = i & 3
		if v == 3 || v == 0 {
			continue
		}
		_, _ = tree.Insert(i, uint8(k), v)
	}

	k = 6
	n = uint64(1 << (k * 2))

	for i = 0; i < n; i++ {
		v = (i >> 4) & 3
		if v == 3 {
			continue
		}
		_, _ = tree.Insert(i, uint8(k), v)
	}

	// tree.Walk(func(code uint64, k uint8, v []uint64) bool {
	// 	t.Logf("%s, %v\n", kmers.Decode(code, int(k)), v)
	// 	return false
	// })

	query := "ACTG"
	t.Logf("query: %s\n", query)
	code, _ := kmers.Encode([]byte(query))
	r, ok = tree.Get(code, uint8(len(query)))
	t.Logf("  %s, %v, %v\n", query, r, ok)

	query = "ACTGC"
	t.Logf("query: %s\n", query)
	code, _ = kmers.Encode([]byte(query))
	_code, _k, _r, ok := tree.LongestPrefix(code, uint8(len(query)))
	if ok {
		t.Logf("  %s, %v, %v\n", kmers.Decode(_code, int(_k)), _r, ok)
	}

	query = "ACGCA"
	code, _ = kmers.Encode([]byte(query))
	srs, _ := tree.Search(code, uint8(len(query)), 4)
	t.Logf("query: %s\n", query)
	for _, sr := range srs {
		t.Logf("  %s, len(prefix): %d, %v\n",
			kmers.Decode(sr.Kmer, int(sr.K)), sr.LenPrefix, sr.Values)
	}

	query = "ACGT"
	code, _ = kmers.Encode([]byte(query))
	nodes, bases := tree.Path(code, uint8(len(query)))
	t.Logf("path of %s: %s, visited nodes: %d, matched bases: %d\n", query, strings.Join(nodes, "->"), len(nodes), bases)
}
