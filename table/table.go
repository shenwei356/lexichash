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

package table

import (
	"encoding/binary"
	"errors"

	"github.com/shenwei356/lexichash/tree"
)

type L3Table struct {
	// K value
	K uint8

	// bps of suffixes as key in L1
	lenS uint8

	// row: 11 bytes
	//	6 bytes:
	//		the last 2 bits: the FLAG of multiple k-mers
	//		content of leading 46 bits depends on the value of FLAG:
	//			flag 0: left bases (k-mer with the suffix removed)
	//			flag 1: 24bits bytes for offset in L2, 22 bytes for number of k-mers
	//	2 bytes: the number of records in V
	//	3 bytes: offset in V
	//
	// nrows: 1<<(2*lenSuffix)
	L1 []byte

	// row: 11 bytes
	//	6 bytes: left bases (k-mer with the suffix removed)
	//	2 bytes: the number of records in V
	//	3 bytes: offset in V
	// nrows: variable, the bigger the LenS, the smaller
	L2 []byte

	// row: 8 bytes, so it would also be uint64
	//	26 bits: ref idx
	//	36 bits: k-mer position
	// 	2 bits: strand flag, 1 for negative strand
	V []uint64
}

var be = binary.BigEndian

// ErrKOverflow means K > 31.
var ErrKOverflow = errors.New("l3table: k-mer size (4-31) overflow, recommended value: [8,10]")

// FromTree creates a L3Table from a k-mer radix tree.
// We assume that all k-mers in the tree have the same value of k
// A lenSuffix value in the range of [8, 10] is recommended.
func FromTree(tree *tree.Tree, lenSuffix uint8) (*L3Table, error) {
	if lenSuffix >= 31 {
		return nil, ErrKOverflow
	}

	var suffixMask uint64 = 1<<(lenSuffix*2) - 1

	// suffix -> kmer -> []uint64
	data := make(map[uint64]map[uint64][]uint64)
	var m map[uint64][]uint64
	first := true
	var K uint8
	var suffix uint64
	var ok bool
	var nrowV int
	tree.Walk(func(code uint64, k uint8, v []uint64) bool {
		if first {
			K = k
			first = false
		}
		suffix = code & suffixMask
		if m, ok = data[suffix]; !ok {
			data[suffix] = map[uint64][]uint64{code: v}
		} else {
			m[code] = v
		}

		nrowV += len(v)

		return false
	})
	var nrowL2 int
	var n int
	for _, m = range data {
		n = len(m)
		if n > 1 {
			nrowL2 += n
		}
	}

	l3t := &L3Table{K: K, lenS: lenSuffix}
	l3t.L1 = make([]byte, 11*(1<<(lenSuffix*2)))
	l3t.L2 = make([]byte, 0, 11*nrowL2)
	l3t.V = make([]uint64, 0, nrowV)

	return l3t, nil
}
