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

package lexichash

import (
	"errors"
	"math"
	"math/rand"
	"sync"

	"github.com/shenwei356/lexichash/iterator"
)

// ErrKOverflow means K > 32.
var ErrKOverflow = errors.New("lexichash: k-mer size overflow, valid range is [4-32]")

// ErrInsufficientMasks means the number of masks is too small.
var ErrInsufficientMasks = errors.New("lexichash: insufficient masks (should be >=4)")

// LexicHash is for finding shared substrings between nucleotide sequences.
type LexicHash struct {
	K int // max length of shared substrings, should be in range of [4, 31]

	Seed  int64    // seed for generating masks
	Masks []uint64 // masks/k-mers

	// pools for storing Kmers and Locs in Mask().
	// users need to call RecycleMaskResult() after using them.
	poolKmers  *sync.Pool
	poolLocses *sync.Pool

	// a []uint64 for storing hashes in Mask(),
	// the object is recycled in Mask().
	poolHashes *sync.Pool
}

// New returns a new LexicHash object.
// nMasks >= 1000 is recommended.
func New(k int, nMasks int) (*LexicHash, error) {
	return NewWithSeed(k, nMasks, 1)
}

// NewWithSeed creates a new LexicHash object with given seed.
// nMasks >= 1000 is recommended.
func NewWithSeed(k int, nMasks int, seed int64) (*LexicHash, error) {
	if k < 4 || k > 32 {
		return nil, ErrKOverflow
	}
	if nMasks < 4 {
		return nil, ErrInsufficientMasks
	}

	lh := &LexicHash{K: k, Seed: seed}

	// ------------ generate masks ------------

	// generate more masks in case there are some duplicated masks
	n := int(float64(nMasks) * 1.2)
	masks := make([]uint64, n)

	r := rand.New(rand.NewSource(seed))
	shift := 64 - k*2
	var mask uint64
	for i := range masks {
		mask = hash64(r.Uint64()) >> shift // hash a random int and cut into k*2 bits
		masks[i] = mask
	}

	uniqUint64s(&masks) // remove duplicates
	if len(masks) > nMasks {
		masks = masks[:nMasks]
	}

	lh.Masks = masks

	// ------------ pools ------------

	lh.poolKmers = &sync.Pool{New: func() interface{} {
		kmers := make([]uint64, len(masks))
		return &kmers
	}}
	lh.poolLocses = &sync.Pool{New: func() interface{} {
		locses := make([][]int, len(masks))
		for i := range locses {
			locses[i] = make([]int, 1)
		}
		return &locses
	}}
	lh.poolHashes = &sync.Pool{New: func() interface{} {
		hashes := make([]uint64, len(masks))
		return &hashes
	}}

	return lh, nil
}

// RecycleMaskResult recycles the results of Mask().
// Please do not forget to call this method after using the mask results.
func (lh *LexicHash) RecycleMaskResult(kmers *[]uint64, locses *[][]int) {
	if kmers != nil {
		lh.poolKmers.Put(kmers)
	}
	if locses != nil {
		lh.poolLocses.Put(locses)
	}
}

// Mask computes the most similar substrings for each mask in sequence s.
// It returns
//
//  1. the list of the most similar k-mers for each mask.
//  2. the start positions of all k-mers, with the last 2 bits as the strand
//     flag (1 for negative strand).
func (lh *LexicHash) Mask(s []byte) (*[]uint64, *[][]int, error) {
	// the k-mer iterator is different from that in
	// https://github.com/shenwei356/bio/blob/master/sketches/iterator.go
	// this one only supports k<=31, with the last two bits as a flag
	// for marking if the k-mer is from the negative strand.
	iter, err := iterator.NewKmerIterator(s, lh.K)
	if err != nil {
		return nil, nil, err
	}

	var masks = lh.Masks
	_kmers := lh.poolKmers.Get().(*[]uint64)  // matched k-mers
	locses := lh.poolLocses.Get().(*[][]int)  // locations of the matched k-mers
	hashes := lh.poolHashes.Get().(*[]uint64) // hashes of matched k-mers
	for i := range *hashes {
		(*hashes)[i] = math.MaxUint64
	}
	var mask, hash, hashRC, h uint64
	var kmer, kmerRC uint64
	var ok bool
	var i, j, js int
	var locs *[]int

	for {
		kmer, kmerRC, ok, _ = iter.NextKmer()
		if !ok {
			break
		}
		j = iter.Index()
		js = j << 2

		for i, mask = range masks {
			h = (*hashes)[i]

			hash = kmer ^ mask
			hashRC = kmerRC ^ mask
			if hashRC < hash { // choose the negative strand
				hash = hashRC
				kmer = kmerRC
				js |= 1 // add the strand flag to the location
			}

			if hash > h {
				continue
			}

			// hash <= h
			locs = &(*locses)[i]
			if hash < h {
				*locs = (*locs)[:1]
				(*locs)[0] = js

				(*hashes)[i] = hash
				(*_kmers)[i] = kmer
			} else {
				*locs = append(*locs, js)
			}
		}
	}

	lh.poolHashes.Put(hashes)
	return _kmers, locses, nil
}
