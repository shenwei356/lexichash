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
	"math/rand"
	"strings"
	"testing"
	"unsafe"

	"github.com/shenwei356/kmers"
	"github.com/shenwei356/lexichash"
)

func TestStructSize(t *testing.T) {
	t.Logf("struct: Sizeof, Alignof\n")
	t.Logf("leafNode: %d, %d", unsafe.Sizeof(leafNode{}), unsafe.Alignof(leafNode{}))
	t.Logf("node: %d, %d", unsafe.Sizeof(node{}), unsafe.Alignof(node{}))
}

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
		tree.Insert(i, uint8(k), v)
	}

	k = 6
	n = uint64(1 << (k * 2))

	for i = 0; i < n; i++ {
		v = (i >> 4) & 3
		if v == 3 {
			continue
		}
		tree.Insert(i, uint8(k), v)
	}

	// t.Logf("number of edges: %d\n", tree.NumEdges())
	t.Logf("number of nodes: %d\n", tree.NumNodes())
	t.Logf("number of leaf nodes: %d\n", tree.NumLeafNodes())

	// tree.Walk(func(code uint64, k uint8, v []uint64) bool {
	// 	t.Logf("%s, %v\n", lexichash.MustDecode(code, k), v)
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
		t.Logf("  %s, %v, %v\n", lexichash.MustDecode(_code, _k), _r, ok)
	}

	query = "ACGCA"
	code, _ = kmers.Encode([]byte(query))
	srs, _ := tree.Search(code, uint8(len(query)), 4)
	t.Logf("query: %s\n", query)
	for _, sr := range *srs {
		t.Logf("  %s, len(prefix): %d, %v\n",
			lexichash.MustDecode(sr.Kmer, sr.K), sr.LenPrefix, sr.Values)
	}
	tree.RecycleSearchResult(srs)

	query = "ACGT"
	code, _ = kmers.Encode([]byte(query))
	nodes, bases := tree.Path(code, uint8(len(query)), uint8(len(query)))
	t.Logf("path of %s: %s, visited nodes: %d, matched bases: %d\n", query, strings.Join(nodes, "->"), len(nodes), bases)
}

func TestBigTree(t *testing.T) {
	nMasks := 1000

	trees := make([]*Tree, nMasks)

	for i := 0; i < nMasks; i++ {
		tree := New()

		k := 21
		nKmers := 10000 // or the number of sequences
		var seed int64 = 1
		seed = int64(i)

		n := int(float64(nKmers) * 1.2)
		codes := make([]uint64, n)

		r := rand.New(rand.NewSource(seed))
		shift := 64 - k*2
		var kmer uint64
		for i := range codes {
			kmer = hash64(r.Uint64()) >> shift // hash a random int and cut into k*2 bits
			codes[i] = kmer
		}
		uniqUint64s(&codes) // remove duplicates
		if len(codes) > nKmers {
			codes = codes[:nKmers]
		}

		for i, kmer := range codes {
			tree.Insert(kmer, uint8(k), uint64(i))
		}

		trees[i] = tree

	}

	tree := trees[0]

	// tree.Walk(func(code uint64, k uint8, v []uint64) bool {
	// 	fmt.Printf("%s\n",lexichash.MustDecode(code, k))
	// 	return false
	// })

	// t.Logf("number of edges: %d\n", tree.NumEdges())
	t.Logf("number of nodes: %d\n", tree.NumNodes())
	t.Logf("number of leaf nodes: %d\n", tree.NumLeafNodes())

	query := "TCATTTAAGTCCAATCGTCAG"
	code, _ := kmers.Encode([]byte(query))
	srs, ok := tree.Search(code, uint8(len(query)), uint8(len(query)))
	t.Logf("query: %s\n", query)
	if ok {
		for _, sr := range *srs {
			t.Logf("  %s, len(prefix): %d, %v\n",
				lexichash.MustDecode(sr.Kmer, sr.K), sr.LenPrefix, sr.Values)
		}
		tree.RecycleSearchResult(srs)
	}

	code, _ = kmers.Encode([]byte(query))
	nodes, bases := tree.Path(code, uint8(len(query)), uint8(len(query)))
	t.Logf("path of %s: %s, visited nodes: %d, matched bases: %d\n", query, strings.Join(nodes, "->"), len(nodes), bases)

}
