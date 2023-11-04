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
	"errors"
	"math"
	"math/bits"
	"math/rand"
	"sort"
	"sync"

	"github.com/shenwei356/kmers"
	iterater "github.com/shenwei356/lexichash/kmer-iterator"
)

// ErrKOverflow means K >= 32.
var ErrKOverflow = errors.New("lexichash: k-mer size (4-31) overflow")

// ErrInsufficientMasks means the number of masks is too small
var ErrInsufficientMasks = errors.New("lexichash: insufficient masks (>=4)")

// LexicHash is for finding shared substrings between nucleotide sequences.
type LexicHash struct {
	K int // max lenth of shared substrings, should be in range of [4, 31]

	Seed  int64    // seed for generating masks
	Masks []uint64 // masks

	// pool for storing Kmers and Locs in Mask().
	// users need to call
	poolKmers *sync.Pool
	poolLocs  *sync.Pool

	// a []uint64 for storing hashes in Mask(),
	// the object is recycled in Mask()
	poolHashes *sync.Pool
	canonical  bool
}

// New returns a new LexicHash object.
// nMasks >= 1000 is recommended.
// Setting canonicalKmer to true is recommended,
// cause it would produces more results.
func New(k int, nMasks int, canonicalKmer bool) (*LexicHash, error) {
	return NewWithSeed(k, nMasks, canonicalKmer, 1)
}

// NewWithSeed creates a new LexicHash object with given seed.
// nMasks >= 1000 is recommended.
// Setting canonicalKmer to true is recommended,
// cause it would produces more results.
func NewWithSeed(k int, nMasks int, canonicalKmer bool, seed int64) (*LexicHash, error) {
	if k < 4 || k > 31 {
		return nil, ErrKOverflow
	}
	if nMasks < 4 {
		return nil, ErrInsufficientMasks
	}

	lh := &LexicHash{K: k, Seed: seed, canonical: canonicalKmer}

	// ------------ generate masks ------------

	// generate more masks in case there are some duplicated masks
	n := int(float64(nMasks) * 1.2)
	masks := make([]uint64, n)

	r := rand.New(rand.NewSource(seed))
	shift := 64 - k*2
	var mask, mask_rc uint64
	for i := range masks {
		mask = hash64(r.Uint64()) >> shift // hash a random int and cut into k*2 bits
		masks[i] = mask

		if lh.canonical {
			mask_rc = kmers.MustRevComp(mask, k) // reverse complement sequence
			if mask_rc < mask {
				masks[i] = mask_rc
			}
		}
	}

	uniqUint64s(&masks) // remove duplicates
	if len(masks) > nMasks {
		masks = masks[:nMasks]
	}

	lh.Masks = masks

	lh.poolKmers = &sync.Pool{New: func() interface{} {
		kmers := make([]uint64, len(masks))
		return &kmers
	}}
	lh.poolLocs = &sync.Pool{New: func() interface{} {
		locs := make([]int, len(masks))
		return &locs
	}}
	lh.poolHashes = &sync.Pool{New: func() interface{} {
		hashes := make([]uint64, len(masks))
		return &hashes
	}}

	return lh, nil
}

// Compare computes the shared substrings, with identical and nested regions removed.
// Each pair of shared substring is recorded as [4]uint64:
//
//	0: k-mer code of shared substring
//	1: length of shared substring
//	2: (start position in s1)<<2 + (negative strand flag)
//	3: (start position in s2)<<2 + (negative strand flag)
//
// Note that the start position is 0-based.
// And the 1-based location of a substring would be
//
//	sub[2]>>2+1, sub[2]>>2+sub[1]
func (lh *LexicHash) Compare(s1 []byte, s2 []byte, minLen int) ([][4]uint64, error) {
	var kmers1, kmers2 *[]uint64
	var locs1, locs2 *[]int
	var err1, err2 error

	var wg sync.WaitGroup
	wg.Add(2)
	go func() {
		kmers1, locs1, err1 = lh.Mask(s1)
		wg.Done()
	}()
	go func() {
		kmers2, locs2, err2 = lh.Mask(s2)
		wg.Done()
	}()
	wg.Wait()
	defer func() {
		lh.RecycleMaskResult(kmers1, locs1)
		lh.RecycleMaskResult(kmers2, locs2)
	}()

	if err1 != nil {
		return nil, err1
	}
	if err2 != nil {
		return nil, err2
	}

	var k2 uint64
	k := lh.K
	delta := k - 32
	var n int // the number of same leading bases
	var loc1, loc2 int
	subs := make([][4]uint64, 0, 8)
	for i, k1 := range *kmers1 {
		k2 = (*kmers2)[i]
		// compute the longest commom prefix.
		// here they have the same K values, so the code below is fine.
		// or you have to align them first.
		// n = int(tree.KmerLongestPrefix(k1>>2, k2>>2, uint8(k), uint8(k)))
		n = bits.LeadingZeros64((k1>>2)^(k2>>2))/2 + delta
		if n < minLen {
			continue
		}

		loc1 = (*locs1)[i]
		if k1&1 > 0 { // it's from the negative strand
			loc1 = loc1 + k - n
		}
		loc2 = (*locs2)[i]
		if k2&1 > 0 { // it's from the negative strand
			loc2 = loc2 + k - n
		}

		subs = append(subs, [4]uint64{
			k1 >> ((k - n + 1) << 1),       // k-mer code of shared substring
			uint64(n),                      // length of shared substring
			uint64(uint64(loc1)<<2 | k1&1), // (start position in s1)<<2 + (negative strand flag)
			uint64(uint64(loc2)<<2 | k2&1), // (start position in s2)<<2 + (negative strand flag)
		})
	}

	// sort and remove duplicates or nested regions
	sort.Slice(subs, func(i, j int) bool {
		if subs[i][2] == subs[j][2] {
			return subs[i][2]>>2+subs[i][1] > subs[j][2]>>2+subs[j][1]
		}
		return subs[i][2] < subs[j][2]
	})

	var i, j int
	var p, v [4]uint64
	var flag bool
	p = subs[0]
	for i = 1; i < len(subs); i++ {
		v = subs[i]
		if v == p || v[2]>>2+v[1] <= p[2]>>2+p[1] { // the same or nested
			if !flag {
				j = i // mark insertion position
				flag = true
			}
			continue
		}

		if flag { // need to insert to previous position
			subs[j] = v
			j++
		}
		p = v
	}
	if j > 0 {
		subs = subs[:j]
	}

	// for _, sub := range subs {
	// 	fmt.Printf("(%3d,%3d, %c) vs (%3d,%3d, %c)  %3d %s\n",
	// 		sub[2]>>2+1, sub[2]>>2+sub[1], strands[sub[2]&1],
	// 		sub[3]>>2+1, sub[3]>>2+sub[1], strands[sub[3]&1],
	// 		sub[1], Kmer2dna(sub[0], int(sub[1])))
	// }

	return subs, nil
}

// RecycleMaskResult recycles the results of Mask().
// Please do not forget to call this method!
func (lh *LexicHash) RecycleMaskResult(kmers *[]uint64, locs *[]int) {
	if kmers != nil {
		lh.poolKmers.Put(kmers)
	}
	if locs != nil {
		lh.poolLocs.Put(locs)
	}
}

// Mask computes the most similar substrings for each mask in sequence s.
// It returns
//
//  1. the list of the most similar k-mers for each mask, with the last 2 bits a strand
//     flag (1 for negative strand)
//  2. the start positions of all k-mers.
//     if k-mer comes from multiple locations, only the smaller one is kept.
func (lh *LexicHash) Mask(s []byte) (*[]uint64, *[]int, error) {
	// the k-mer iterator is different from that in
	// https://github.com/shenwei356/bio/blob/master/sketches/iterator.go
	// this one only supports k<=31, with the last two bits as a flag
	// for marking if the k-mer is from the negative strand.
	iter, err := iterater.NewKmerIterator(s, lh.K, lh.canonical)
	if err != nil {
		return nil, nil, err
	}

	var masks = lh.Masks
	// _kmers := make([]uint64, len(masks)) // matched k-mers
	// hashes := make([]uint64, len(masks)) // hashes of matched k-mers
	// locs := make([]int, len(masks))     // locations of the matched k-mers
	_kmers := lh.poolKmers.Get().(*[]uint64)
	locs := lh.poolLocs.Get().(*[]int)
	hashes := lh.poolHashes.Get().(*[]uint64) // hashes of matched k-mers
	for i := range *hashes {
		(*hashes)[i] = math.MaxUint64
	}

	var mask, hash uint64
	var kmer [2]uint64
	var ok bool
	var i, j int
	canonical := lh.canonical

	if canonical {
		for {
			kmer, ok, _ = iter.NextKmer()
			if !ok {
				break
			}
			j = iter.Index()

			for i, mask = range masks {
				hash = kmer[0]>>2 ^ mask
				// smaller hash
				if hash < (*hashes)[i] {
					(*hashes)[i] = hash
					(*_kmers)[i] = kmer[0]
					(*locs)[i] = j
				}
			}
		}

		lh.poolHashes.Put(hashes)
		return _kmers, locs, nil
	}

	for {
		kmer, ok, _ = iter.NextKmer()
		if !ok {
			break
		}
		j = iter.Index()

		for i, mask = range masks {
			hash = kmer[0]>>2 ^ mask
			// smaller hash
			if hash < (*hashes)[i] {
				(*hashes)[i] = hash
				(*_kmers)[i] = kmer[0]
				(*locs)[i] = j
			}

			// try both strands
			hash = kmer[1]>>2 ^ mask
			// smaller hash
			if hash < (*hashes)[i] {
				(*hashes)[i] = hash
				(*_kmers)[i] = kmer[1]
				(*locs)[i] = j
			}
		}
	}

	lh.poolHashes.Put(hashes)
	return _kmers, locs, nil
}
