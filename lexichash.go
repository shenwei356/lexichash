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
	"sort"
	"sync"

	"github.com/shenwei356/lexichash/iterator"
)

// ErrKOverflow means K > 32.
var ErrKOverflow = errors.New("lexichash: k-mer size overflow, valid range is [3-32]")

// ErrInsufficientMasks means the number of masks is too small.
var ErrInsufficientMasks = errors.New("lexichash: insufficient masks (should be >=4)")

// LexicHash is for finding shared substrings between nucleotide sequences.
type LexicHash struct {
	K int // max length of shared substrings, should be in range of [4, 31]

	Seed  int64    // seed for generating masks
	Masks []uint64 // masks/k-mers

	// indexes for fast locating masks to compare.
	// m3 means using the first 3 bases as the map keys,
	// Three has the best balance in tests using 5000 and 10000 masks,
	// smaller or bigger values would slow down the speed significantly.
	m3 []*[]int
	// pool for checking masks without matches.
	// sync.Pool is used because Mask() method might be called concurrently.
	// the object is recycled in Mask().
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
// p is the length of mask k-mer prefixes which need to be checked for low-complexity.
// p == 0 for no checking.
func New(k int, nMasks int, p int) (*LexicHash, error) {
	return NewWithSeed(k, nMasks, 1, p)
}

// NewWithSeed creates a new LexicHash object with given seed.
// nMasks better be >= 1024 and better be power of 4,
// i.e., 4, 16, 64, 256, 1024, 4096 ...
// p is the length of mask k-mer prefixes which need to be checked for low-complexity.
// p == 0 for no checking.
func NewWithSeed(k int, nMasks int, randSeed int64, p int) (*LexicHash, error) {
	if k < 3 || k > 32 {
		return nil, ErrKOverflow
	}
	if nMasks < 4 {
		return nil, ErrInsufficientMasks
	}

	lh := &LexicHash{K: k, Seed: randSeed}

	// ------------ generate masks ------------

	masks := genRandomMasks(k, nMasks, randSeed, p)

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

// p is the length of prefixes which need to be checked for low-complexity.
func genRandomMasks(k int, nMasks int, randSeed int64, p int) []uint64 {
	masks := make([]uint64, nMasks)
	m := make(map[uint64]interface{}, nMasks) // to avoid duplicates
	r := rand.New(rand.NewSource(randSeed))
	checkLC := p > 0

	// generate 4^x prefix
	nPrefix := 1
	for 1<<(nPrefix<<1) <= nMasks {
		nPrefix++
	}
	nPrefix--
	n := 1 << (nPrefix << 1)
	var bases []uint64
	if checkLC && nPrefix >= 5 { // filter out k-mer with a prefix of AAAAA, CCCCC, GGGGG or TTTTT
		bases = make([]uint64, 0, n)
		for i := 0; i < n; i++ {
			if IsLowComplexity(uint64(i), nPrefix) {
				continue
			}
			bases = append(bases, uint64(i))
		}
		n = len(bases)
	} else {
		bases = make([]uint64, n)
		for i := 0; i < n; i++ {
			bases[i] = uint64(i)
		}
	}

	// distribute these prefixes
	var j = 0
	for j = 0; j < nMasks/n; j++ {
		copy(masks[j*n:(j+1)*n], bases)
	}
	if nMasks%n != 0 { // randomly sampling for left
		r.Shuffle(n, func(i, j int) { bases[i], bases[j] = bases[j], bases[i] })
		copy(masks[j*n:], bases[:nMasks%n])
	}

	// concatenate with random numbers
	var _mask uint64 = 1<<(uint64(k-nPrefix)<<1) - 1
	shiftP := uint64(k-nPrefix) << 1
	var mask uint64
	var v uint64
	var i int
	var ok bool
	var prefix uint64
	var tries int
	for {
		v = r.Uint64()
		mask = Hash64(v)&_mask | masks[i]<<shiftP

		if checkLC {
			prefix = mask >> ((k - p) << 1)
			if IsLowComplexity(prefix, p) {
				tries++
				if tries <= 10 {
					continue
				}
				// give up trying
			}
		}

		if _, ok = m[mask]; ok {
			continue
		}
		masks[i] = mask
		m[mask] = struct{}{}
		i++
		tries = 0

		if i == nMasks {
			break
		}
	}

	// sort
	sort.Slice(masks, func(i, j int) bool { return masks[i] < masks[j] })

	return masks
}

// indexMasks indexes masks with lists for fast locating masks to compare
func (lh *LexicHash) indexMasks() {
	k := lh.K
	var prefix uint64
	var list *[]int

	// 3
	m := make([]*[]int, 1<<(3<<1))
	for i, mask := range lh.Masks {
		prefix = mask >> ((k - 3) << 1)
		if list = m[prefix]; list == nil {
			m[prefix] = &[]int{i}
		} else {
			*list = append(*list, i)
		}
	}
	lh.m3 = m
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
//
// skipRegions is optional, which is used to skip some masked regions.
// E.g., in reference indexing step, contigs of a genome can be concatenated with k-1 N's,
// where need to be ommitted.
//
// The regions should be 0-based and ascendingly sorted.
// e.g., [100, 130], [200, 230] ...
func (lh *LexicHash) Mask(s []byte, skipRegions [][2]int) (*[]uint64, *[][]int, error) {
	_kmers := lh.poolKmers.Get().(*[]uint64)  // matched k-mers
	locses := lh.poolLocses.Get().(*[][]int)  // locations of the matched k-mers
	hashes := lh.poolHashes.Get().(*[]uint64) // hashes of matched k-mers
	for i := range *hashes {
		(*hashes)[i] = math.MaxUint64
	}

	masks := lh.Masks
	k := lh.K
	var mask, hash, h uint64
	var kmer, kmerRC uint64
	var ok bool
	var i, j, js int
	var locs *[]int

	nRegions := len(skipRegions)
	checkRegion := nRegions > 0
	var ri, rs, re int

	// -----------------------------------------------------------------------------
	// round 1: fast locating of mask to compare with prefix indexing

	// This k-mer iterator is a simplified version of
	// https://github.com/shenwei356/bio/blob/master/sketches/iterator.go
	iter, err := iterator.NewKmerIterator(s, lh.K)
	if err != nil {
		return nil, nil, err
	}

	m3 := lh.m3

	if checkRegion {
		ri = 0
		rs, re = skipRegions[ri][0], skipRegions[ri][1]
	}

	for {
		kmer, kmerRC, ok, _ = iter.NextKmer()
		if !ok {
			break
		}
		if kmer == 0 { // all bases are A's or N's.
			continue
		}

		j = iter.Index()

		// skip some regions
		if checkRegion && rs <= j {
			if j <= re { // in the region
				if j == re { // update region to check
					ri++
					if ri == nRegions { // this is already the last one
						checkRegion = false
					} else {
						rs, re = skipRegions[ri][0], skipRegions[ri][1]
					}
				}

				continue
			}
		}

		js = j << 2

		// ---------- positive strand ----------

		for _, i = range *m3[int(kmer>>((k-3)<<1))] {
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

		for _, i = range *m3[int(kmerRC>>((k-3)<<1))] {
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

	// -----------------------------------------------------------------------------
	// round 2.
	// some masks may not have any matches,
	// use the classic method for them, i.e., compare with all left masks one by one.

	noMatches := lh.poolList.Get().(*[]int)
	*noMatches = (*noMatches)[:0]
	for i, h := range *hashes {
		if h == math.MaxUint64 {
			*noMatches = append(*noMatches, i)
		}
	}

	if len(*noMatches) == 0 { // cool, no need to continue
		lh.poolList.Put(noMatches)
		lh.poolHashes.Put(hashes)
		return _kmers, locses, nil
	}

	// don't worry the efficiency, NewKmerIterator is optimized to reuse objects
	iter, _ = iterator.NewKmerIterator(s, lh.K)

	if checkRegion {
		ri = 0
		rs, re = skipRegions[ri][0], skipRegions[ri][1]
	}

	for {
		kmer, kmerRC, ok, _ = iter.NextKmer()
		if !ok {
			break
		}

		if kmer == 0 { // all bases are A's or N's.
			continue
		}

		j = iter.Index()

		// skip some regions
		if checkRegion && rs <= j {
			if j <= re { // in the region
				if j == re { // update region to check
					ri++
					if ri == nRegions { // this is already the last one
						checkRegion = false
					} else {
						rs, re = skipRegions[ri][0], skipRegions[ri][1]
					}
				}
				continue
			}
		}

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

	// recycle some reusable objects
	lh.poolHashes.Put(hashes)
	lh.poolList.Put(noMatches)
	return _kmers, locses, nil
}
