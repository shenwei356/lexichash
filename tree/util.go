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
	"math/bits"

	"github.com/shenwei356/kmers"
	"github.com/twotwotwo/sorts/sortutil"
)

var bit2base = [4]byte{'A', 'C', 'G', 'T'}

func Kmer2dna(code uint64, k int) []byte {
	return kmers.Decode(code, k)
}

// KmerBaseAt returns the base in pos i (0-based).
func KmerBaseAt(code uint64, k uint8, i uint8) uint8 {
	return uint8(code >> ((k - i - 1) << 1) & 3)
}

// KmerPrefix returns the first n bases. n needs to be > 0.
// The length of the prefix is n.
func KmerPrefix(code uint64, k uint8, n uint8) uint64 {
	return code >> ((k - n) << 1)
}

// KmerSuffix returns the suffix starting from position i (0-based).
// The length of the suffix is k - commonPrefix.
func KmerSuffix(code uint64, k uint8, i uint8) uint64 {
	return code & (1<<((k-i)<<1) - 1)
}

// KmerLongestPrefix returns the length of the longest prefix.
func KmerLongestPrefix(code1, code2 uint64, k1, k2 uint8) uint8 {
	var d uint8
	if k1 >= k2 { // most of the cases
		code1 = code1 >> ((k1 - k2) << 1)
		d = 32 - k2
	} else {
		code2 = code2 >> ((k2 - k1) << 1)
		d = 32 - k1
	}
	return uint8(bits.LeadingZeros64(code1^code2)>>1) - d
}

// KmerHasPrefix check if a k-mer has a prefix
func KmerHasPrefix(code uint64, prefix uint64, k1, k2 uint8) bool {
	if k1 < k2 {
		return false
	}
	return code>>((k1-k2)<<1) == prefix
}

// MustKmerHasPrefix check if a k-mer has a prefix, without if branch.
func MustKmerHasPrefix(code uint64, prefix uint64, k1, k2 uint8) bool {
	return code>>((k1-k2)<<1) == prefix
}

// https://gist.github.com/badboy/6267743 .
// version with mask: https://gist.github.com/lh3/974ced188be2f90422cc .
func hash64(key uint64) uint64 {
	key = (^key) + (key << 21) // key = (key << 21) - key - 1
	key = key ^ (key >> 24)
	key = (key + (key << 3)) + (key << 8) // key * 265
	key = key ^ (key >> 14)
	key = (key + (key << 2)) + (key << 4) // key * 21
	key = key ^ (key >> 28)
	key = key + (key << 31)
	return key
}

func uniqUint64s(list *[]uint64) {
	if len(*list) == 0 || len(*list) == 1 {
		return
	}

	sortutil.Uint64s(*list)

	var i, j int
	var p, v uint64
	var flag bool
	p = (*list)[0]
	for i = 1; i < len(*list); i++ {
		v = (*list)[i]
		if v == p {
			if !flag {
				j = i // mark insertion position
				flag = true
			}
			continue
		}

		if flag { // need to insert to previous position
			(*list)[j] = v
			j++
		}
		p = v
	}
	if j > 0 {
		*list = (*list)[:j]
	}
}
