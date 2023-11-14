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

package iterater

import (
	"errors"
	"fmt"
	"sync"

	"github.com/shenwei356/kmers"
)

// ErrInvalidK means k < 1.
var ErrInvalidK = fmt.Errorf("k-mer iterator: invalid k-mer size (1 <= k <= 31)")

// ErrEmptySeq sequence is empty.
var ErrEmptySeq = fmt.Errorf("k-mer iterator: empty sequence")

// ErrShortSeq means the sequence is shorter than k
var ErrShortSeq = fmt.Errorf("k-mer iterator: sequence too short")

// ErrIllegalBase means that base beyond IUPAC symbols are  detected.
var ErrIllegalBase = errors.New("k-mer iterator: illegal base")

// ErrKTooLarge means that the k-mer size is too large.
var ErrKTooLarge = fmt.Errorf("k-mer iterator: k-mer size is too large")

// ErrInvalidM means that the m-mer size is too large or too small, should be in range of [4, k].
var ErrInvalidM = fmt.Errorf("k-mer iterator: invalid m-mer size, should be in range of [4, k]")

// ErrInvalidScale means
var ErrInvalidScale = fmt.Errorf("k-mer iterator: invalid scale, should be in range of [1, k-m+1]")

var poolIterator = &sync.Pool{New: func() interface{} {
	return &Iterator{}
}}

// Iterator is a kmer code (k<=31) or hash iterator.
type Iterator struct {
	s       []byte
	k       int
	kUint   uint // uint(k)
	kP1     int  // k -1
	kP1Uint uint // uint(k-1)

	finished bool
	idx      int

	// for KmerIterator
	canonical bool
	length    int
	end, e    int
	first     bool
	kmer      []byte
	codeBase  uint64
	preCode   uint64
	preCodeRC uint64
	codeRC    uint64

	mask1 uint64 // (1<<(kP1Uint*2))-1
	mask2 uint   // iter.kP1Uint*2
}

// NewKmerIterator returns k-mer code iterator.
func NewKmerIterator(s []byte, k int, canonical bool) (*Iterator, error) {
	if k < 1 {
		return nil, ErrInvalidK
	}
	if len(s) < k {
		return nil, ErrShortSeq
	}

	// iter := &Iterator{s: s2, k: k, circular: circular}
	iter := poolIterator.Get().(*Iterator)
	iter.s = s
	iter.k = k
	iter.finished = false
	iter.idx = 0

	iter.canonical = canonical
	iter.length = len(s)
	iter.end = iter.length - k + 1
	iter.kUint = uint(k)
	iter.kP1 = k - 1
	iter.kP1Uint = uint(k - 1)
	iter.mask1 = (1 << (iter.kP1Uint << 1)) - 1
	iter.mask2 = iter.kP1Uint << 1

	iter.first = true

	return iter, nil
}

// NextKmer returns next two k-mer codes.
// codes[0] is for the positive strand,
// codes[1] is for the negative strand.
func (iter *Iterator) NextKmer() (code, codeRC uint64, ok bool, err error) {
	if iter.finished {
		return 0, 0, false, nil
	}

	if iter.idx == iter.end {
		iter.finished = true
		poolIterator.Put(iter)
		return 0, 0, false, nil
	}

	iter.e = iter.idx + iter.k
	iter.kmer = iter.s[iter.idx:iter.e]

	if !iter.first {
		iter.codeBase = base2bit[iter.kmer[iter.kP1]]
		if iter.codeBase == 4 {
			err = ErrIllegalBase
		}

		// compute code from previous one
		code = (iter.preCode&iter.mask1)<<2 | iter.codeBase

		// compute code of revcomp kmer from previous one
		iter.codeRC = (iter.codeBase^3)<<(iter.mask2) | (iter.preCodeRC >> 2)
	} else {
		code, err = kmers.Encode(iter.kmer)
		iter.codeRC = kmers.MustRevComp(code, iter.k)
		iter.first = false
	}
	if err != nil {
		return 0, 0, false, fmt.Errorf("encode %s: %s", iter.kmer, err)
	}

	iter.preCode = code
	iter.preCodeRC = iter.codeRC
	iter.idx++

	code = code << 2            // positive strand
	codeRC = iter.codeRC<<2 | 1 // negative strand

	if iter.canonical && code > iter.codeRC {
		code = codeRC
	}

	return code, codeRC, true, nil
}

// Index returns current 0-baesd index.
func (iter *Iterator) Index() int {
	return iter.idx - 1
}

var base2bit = [256]uint64{
	4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
	4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
	4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
	4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
	4, 0, 1, 1, 0, 4, 4, 2, 0, 4, 4, 2, 4, 0, 0, 4,
	4, 4, 0, 1, 3, 3, 0, 0, 4, 1, 4, 4, 4, 4, 4, 4,
	4, 0, 1, 1, 0, 4, 4, 2, 0, 4, 4, 2, 4, 0, 0, 4,
	4, 4, 0, 1, 3, 3, 0, 0, 4, 1, 4, 4, 4, 4, 4, 4,
	4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
	4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
	4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
	4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
	4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
	4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
	4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
	4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
}

var bit2base = [4]byte{'A', 'C', 'G', 'T'}
