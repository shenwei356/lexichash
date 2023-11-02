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
	"fmt"
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
		if v == 3 {
			continue
		}
		_, _ = tree.Insert(i, uint8(k), v)
	}
	shortKmer, err := kmers.Encode([]byte("ACT"))
	if err != nil {
		t.Error(err)
		return
	}
	_, _ = tree.Insert(shortKmer, uint8(3), v)

	// tree.Walk(func(code uint64, k uint8, v []uint64) bool {
	// 	fmt.Printf("%s, %v\n", kmers.Decode(code, int(k)), v)
	// 	return false
	// })

	query := "ACTG"
	fmt.Printf("\nquery: %s\n", query)
	code, _ := kmers.Encode([]byte(query))
	r, ok = tree.Get(code, uint8(len(query)))
	fmt.Printf("  %s, %v, %v\n", query, r, ok)

	query = "ACTC"
	fmt.Printf("\nquery: %s\n", query)
	code, _ = kmers.Encode([]byte(query))
	_code, _k, _r, ok := tree.LongestPrefix(code, uint8(len(query)))
	fmt.Printf("  %s, %v, %v\n", kmers.Decode(_code, int(_k)), _r, ok)

	query = "ACGT"
	code, _ = kmers.Encode([]byte(query))
	_codes, _ks, _rs, ok := tree.Search(code, uint8(len(query)), 3)
	fmt.Printf("\nquery: %s\n", query)
	for i, _code := range _codes {
		_k := _ks[i]
		_r := _rs[i]

		fmt.Printf("%s, %v, %v\n", kmers.Decode(_code, int(_k)), _r, ok)
	}
}
