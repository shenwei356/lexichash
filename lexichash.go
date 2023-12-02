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
	"github.com/shenwei356/lexichash/util"
	"github.com/twotwotwo/sorts/sortutil"
)

// ErrKOverflow means K > 32.
var ErrKOverflow = errors.New("lexichash: k-mer size overflow, valid range is [5-32]")

// ErrInsufficientMasks means the number of masks is too small.
var ErrInsufficientMasks = errors.New("lexichash: insufficient masks (should be >=4)")

// LexicHash is for finding shared substrings between nucleotide sequences.
type LexicHash struct {
	K int // max length of shared substrings, should be in range of [4, 31]

	Seed  int64    // seed for generating masks
	Masks []uint64 // masks/k-mers

	// indexes for fast locating masks to compare
	m1 []*[]int
	m2 []*[]int
	m3 []*[]int
	m4 []*[]int
	m5 []*[]int
	// pool for checking masks without matches
	poolList *sync.Pool

	// pools for storing Kmers and Locs in Mask().
	// users need to call RecycleMaskResult() after using them.
	poolKmers  *sync.Pool
	poolLocses *sync.Pool

	// a []uint64 for storing hashes in Mask(),
	// the object is recycled in Mask().
	poolHashes *sync.Pool
}

// New returns a new LexicHash object.
// nMasks better be >= 1024 and better be power of 4,
// i.e., 4, 16, 64, 256, 1024, 4096 ...
func New(k int, nMasks int) (*LexicHash, error) {
	return NewWithSeed(k, nMasks, 1)
}

// NewWithSeed creates a new LexicHash object with given seed.
// nMasks better be >= 1024 and better be power of 4,
// i.e., 4, 16, 64, 256, 1024, 4096 ...
func NewWithSeed(k int, nMasks int, randSeed int64) (*LexicHash, error) {
	if k < 5 || k > 32 {
		return nil, ErrKOverflow
	}
	if nMasks < 4 {
		return nil, ErrInsufficientMasks
	}

	lh := &LexicHash{K: k, Seed: randSeed}

	// ------------ generate masks ------------

	masks := genRandomMasks(k, nMasks, randSeed)

	lh.Masks = masks
	lh.indexMasks()

	// ------------ pools ------------

	lh.poolList = &sync.Pool{New: func() interface{} {
		tmp := make([]int, 128)
		return &tmp
	}}
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

func genRandomMasks(k int, nMasks int, randSeed int64) []uint64 {
	masks := make([]uint64, nMasks)
	m := make(map[uint64]interface{}, nMasks) // to avoid duplicates
	r := rand.New(rand.NewSource(randSeed))
	// var _mask uint64 = 1<<(k<<1) - 1

	// generate 4^x prefix
	nPrefix := 1
	for 1<<(nPrefix<<1) <= nMasks {
		nPrefix++
	}
	nPrefix--
	n := 1 << (nPrefix << 1)
	bases := make([]uint64, n)
	for i := 0; i < n; i++ {
		bases[i] = uint64(i)
	}

	// distribute these prefixes
	var j = 0
	for j = 0; j < nMasks/n; j++ {
		copy(masks[j*n:(j+1)*n], bases)
	}
	if nMasks%n != 0 { // randomly sampling for left
		rand.Shuffle(n, func(i, j int) { bases[i], bases[j] = bases[j], bases[i] })
		copy(masks[j*n:], bases[:nMasks%n])
	}

	// concatenate with random numbers
	var _mask uint64 = 1<<(uint64(k-nPrefix)<<1) - 1
	shiftP := uint64(k-nPrefix) << 1
	var mask uint64
	var v uint64
	var i int
	var ok bool
	for {
		v = r.Uint64()
		mask = util.Hash64(v)&_mask | masks[i]<<shiftP
		if _, ok = m[mask]; ok {
			continue
		}
		masks[i] = mask
		m[mask] = struct{}{}
		i++

		if i == nMasks {
			break
		}
	}

	// sort
	sortutil.Uint64s(masks)

	return masks
}

// indexMasks indexes masks with lists for fast locating masks to compare
func (lh *LexicHash) indexMasks() {
	k := lh.K
	var prefix uint64
	var list *[]int
	// 5
	m := make([]*[]int, 1<<(5<<1))
	for i, mask := range lh.Masks {
		prefix = mask >> ((k - 5) << 1)
		if list = m[prefix]; list == nil {
			m[prefix] = &[]int{i}
		} else {
			*list = append(*list, i)
		}
	}
	lh.m5 = m

	// 4
	m = make([]*[]int, 1<<(4<<1))
	for i, mask := range lh.Masks {
		prefix = mask >> ((k - 4) << 1)
		if list = m[prefix]; list == nil {
			m[prefix] = &[]int{i}
		} else {
			*list = append(*list, i)
		}
	}
	lh.m4 = m

	// 3
	m = make([]*[]int, 1<<(3<<1))
	for i, mask := range lh.Masks {
		prefix = mask >> ((k - 3) << 1)
		if list = m[prefix]; list == nil {
			m[prefix] = &[]int{i}
		} else {
			*list = append(*list, i)
		}
	}
	lh.m3 = m

	// 2
	m = make([]*[]int, 1<<(2<<1))
	for i, mask := range lh.Masks {
		prefix = mask >> ((k - 2) << 1)
		if list = m[prefix]; list == nil {
			m[prefix] = &[]int{i}
		} else {
			*list = append(*list, i)
		}
	}
	lh.m2 = m

	// 1
	m = make([]*[]int, 1<<(1<<1))
	for i, mask := range lh.Masks {
		prefix = mask >> ((k - 1) << 1)
		if list = m[prefix]; list == nil {
			m[prefix] = &[]int{i}
		} else {
			*list = append(*list, i)
		}
	}
	lh.m1 = m
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
//  2. the start 0-based positions of all k-mers, with the last 2 bits as the strand
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

	masks := lh.Masks
	_kmers := lh.poolKmers.Get().(*[]uint64)  // matched k-mers
	locses := lh.poolLocses.Get().(*[][]int)  // locations of the matched k-mers
	hashes := lh.poolHashes.Get().(*[]uint64) // hashes of matched k-mers
	for i := range *hashes {
		(*hashes)[i] = math.MaxUint64
	}
	var mask, hash, h uint64
	var kmer, kmerRC uint64
	var ok bool
	var i, j, js int
	var locs *[]int

	k := lh.K
	m5 := lh.m5
	m4 := lh.m4
	m3 := lh.m3
	m2 := lh.m2
	m1 := lh.m1
	var maskIdxs *[]int

	for {
		kmer, kmerRC, ok, _ = iter.NextKmer()
		if !ok {
			break
		}
		j = iter.Index()
		js = j << 2

		// ---------- positive strand ----------

		maskIdxs = m5[int(kmer>>((k-5)<<1))]
		if maskIdxs == nil {
			maskIdxs = m4[int(kmer>>((k-4)<<1))]
			if maskIdxs == nil {
				maskIdxs = m3[int(kmer>>((k-3)<<1))]
				if maskIdxs == nil {
					maskIdxs = m2[int(kmer>>((k-2)<<1))]
					if maskIdxs == nil {
						maskIdxs = m1[int(kmer>>((k-1)<<1))]
					}
				}
			}
		}

		for _, i = range *maskIdxs {
			mask = masks[i]
			h = (*hashes)[i]

			hash = kmer ^ mask

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

		// ---------- negative strand ----------

		js |= 1 // add the strand flag to the location

		maskIdxs = m5[int(kmerRC>>((k-5)<<1))]
		if maskIdxs == nil {
			maskIdxs = m4[int(kmerRC>>((k-4)<<1))]
			if maskIdxs == nil {
				maskIdxs = m3[int(kmerRC>>((k-3)<<1))]
				if maskIdxs == nil {
					maskIdxs = m2[int(kmerRC>>((k-2)<<1))]
					if maskIdxs == nil {
						maskIdxs = m1[int(kmerRC>>((k-1)<<1))]
					}
				}
			}
		}

		for _, i = range *maskIdxs {
			mask = masks[i]
			h = (*hashes)[i]

			hash = kmerRC ^ mask

			if hash > h {
				continue
			}

			// hash <= h
			locs = &(*locses)[i]
			if hash < h {
				*locs = (*locs)[:1]
				(*locs)[0] = js

				(*hashes)[i] = hash
				(*_kmers)[i] = kmerRC
			} else {
				*locs = append(*locs, js)
			}
		}
	}

	// -------------------------------------------------
	// some mask may do not have any matches,
	// use the classic method for them

	noMatches := lh.poolList.Get().(*[]int)
	*noMatches = (*noMatches)[:0]
	for i, h := range *hashes {
		if h == math.MaxUint64 {
			*noMatches = append(*noMatches, i)
		}
	}

	if len(*noMatches) == 0 { // cool, no need to continue
		lh.poolList.Put(noMatches)
		return _kmers, locses, nil
	}

	iter, _ = iterator.NewKmerIterator(s, lh.K)
	for {
		kmer, kmerRC, ok, _ = iter.NextKmer()
		if !ok {
			break
		}
		j = iter.Index()
		js = j << 2

		// ---------- positive strand ----------

		for _, i = range *noMatches {
			mask = masks[i]
			h = (*hashes)[i]

			hash = kmer ^ mask

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

		// ---------- negative strand ----------

		js |= 1 // add the strand flag to the location

		for _, i = range *noMatches {
			mask = masks[i]
			h = (*hashes)[i]

			hash = kmerRC ^ mask

			if hash > h {
				continue
			}

			// hash <= h
			locs = &(*locses)[i]
			if hash < h {
				*locs = (*locs)[:1]
				(*locs)[0] = js

				(*hashes)[i] = hash
				(*_kmers)[i] = kmerRC
			} else {
				*locs = append(*locs, js)
			}
		}
	}

	lh.poolList.Put(noMatches)
	return _kmers, locses, nil
}
