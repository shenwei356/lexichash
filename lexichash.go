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
	"bufio"
	"bytes"
	"errors"
	"fmt"
	"math"
	"math/rand"
	"os"
	"sort"
	"sync"

	"github.com/shenwei356/kmers"
	"github.com/shenwei356/lexichash/iterator"
)

// ErrKOverflow means K > 32.
var ErrKOverflow = errors.New("lexichash: k-mer size overflow, valid range is [3-32]")

// ErrPrefixOverflow means prefix > k.
var ErrPrefixOverflow = errors.New("lexichash: prefix should be in range of [3, k]")

// ErrInsufficientMasks means the number of masks is too small.
var ErrInsufficientMasks = errors.New("lexichash: insufficient masks (should be >=64)")

// LexicHash is for finding shared substrings between nucleotide sequences.
type LexicHash struct {
	K int // max length of shared substrings, should be in range of [4, 31]

	Seed  int64    // seed for generating masks
	Masks []uint64 // masks/k-mers

	// indexes for fast locating masks to compare.
	// m3 means using the first 3 bases as the map keys.
	// for long sequence, e.g., > 1Mb, any 5-bp prefix probably has some matches.
	// while for short sequences, we use a 3-bp prefix for fast locating.
	// For masks > 64, 3-bp prefix must have a match.
	m3 []*[]int
	m5 []*[]int

	// mN     map[uint64]*[]int // the length of the prefix is given by IndexMasks
	mN     []*[]int // slice is faster than map
	prefix int      //

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
// nMasks should be >=64, and better be >= 1024 and better be power of 4,
// i.e., 64, 256, 1024, 4096 ...
// p is the length of mask k-mer prefixes which need to be checked for low-complexity.
// p == 0 for no checking.
func New(k int, nMasks int, p int) (*LexicHash, error) {
	return NewWithSeed(k, nMasks, 1, p)
}

// NewWithSeed creates a new LexicHash object with given seed.
// nMasks should be >=64, and better be >= 1024 and better be power of 4,
// i.e., 64, 256, 1024, 4096 ...
// p is the length of mask k-mer prefixes which need to be checked for low-complexity.
// p == 0 for no checking.
func NewWithSeed(k int, nMasks int, randSeed int64, p int) (*LexicHash, error) {
	if k < 3 || k > 32 {
		return nil, ErrKOverflow
	}
	if nMasks < 64 { // 4*4*4 = 64
		return nil, ErrInsufficientMasks
	}

	masks := genRandomMasks(k, nMasks, randSeed, p)

	lh, err := NewWithMasks(k, masks)
	if err != nil {
		return lh, err
	}

	lh.Seed = randSeed
	return lh, nil
}

// NewWithMasks creates a new LexicHash object with custom kmers.
// nMasks should be >=64, and better be >= 1024 and better be power of 4,
// i.e., 64, 256, 1024, 4096 ...
func NewWithMasks(k int, masks []uint64) (*LexicHash, error) {
	if k < 3 || k > 32 {
		return nil, ErrKOverflow
	}
	if len(masks) < 64 { // 4*4*4 = 64
		return nil, ErrInsufficientMasks
	}
	// checking mask
	maxMask := 1<<(k<<1) - 1
	for mask := range masks {
		if mask > maxMask {
			return nil, fmt.Errorf("lexichash: given invalid k-mer code for k=%d: %d", k, mask)
		}
	}

	lh := &LexicHash{K: k}

	// sort
	sort.Slice(masks, func(i, j int) bool { return masks[i] < masks[j] })

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

// NewFromTextFile creates a new LexicHash object with custom kmers in a txt file.
func NewFromTextFile(file string) (*LexicHash, error) {
	fh, err := os.Open(file)
	if err != nil {
		return nil, err
	}

	masks := make([]uint64, 0, 1024)

	scanner := bufio.NewScanner(fh)
	var line []byte
	var kmer uint64
	var k int
	for scanner.Scan() {
		line = bytes.Trim(scanner.Bytes(), "\r\n")
		if len(line) == 0 || line[0] == '#' {
			continue
		}
		if k == 0 {
			k = len(line)
		} else if k != len(line) {
			return nil, fmt.Errorf("lexichash: inconsistent k-mer length: %d != %d, %s", k, len(line), line)
		}

		kmer, err = kmers.Encode(line)
		if err != nil {
			return nil, err
		}

		masks = append(masks, kmer)
	}
	if err = scanner.Err(); err != nil {
		return nil, err
	}
	fh.Close()

	return NewWithMasks(k, masks)
}

// p is the length of prefixes which need to be checked for low-complexity.
func genRandomMasks(k int, nMasks int, randSeed int64, p int) []uint64 {
	masks := make([]uint64, nMasks)
	m := make(map[uint64]interface{}, nMasks) // to avoid duplicates
	r := rand.New(rand.NewSource(randSeed))
	if p > k {
		p = k
	}
	checkLC := p > 0

	// generate 4^x prefix
	lenPrefix := 1
	for 1<<(lenPrefix<<1) <= nMasks {
		lenPrefix++
	}
	lenPrefix--
	n := 1 << (lenPrefix << 1)
	var bases []uint64
	if checkLC && lenPrefix >= 5 { // filter out k-mer with a prefix of AAAAA, CCCCC, GGGGG or TTTTT
		bases = make([]uint64, 0, n)
		for i := 0; i < n; i++ {
			if IsLowComplexity(uint64(i), lenPrefix) {
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
	var _mask uint64 = 1<<(uint64(k-lenPrefix)<<1) - 1
	shiftP := uint64(k-lenPrefix) << 1
	var mask uint64
	var v uint64
	var i int
	var ok bool
	var prefix uint64
	var tries int
	shiftOffset := (k - p) << 1
	for {
		v = r.Uint64()
		mask = Hash64(v)&_mask | masks[i]<<shiftP

		if checkLC {
			prefix = mask >> shiftOffset
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

	shiftOffset := (k - 3) << 1
	m := make([]*[]int, 1<<(3<<1))
	for i, mask := range lh.Masks {
		prefix = mask >> shiftOffset
		if list = m[prefix]; list == nil {
			m[prefix] = &[]int{i}
		} else {
			*list = append(*list, i)
		}
	}
	lh.m3 = m

	// 5
	shiftOffset = (k - 5) << 1
	m = make([]*[]int, 1<<(5<<1))
	for i, mask := range lh.Masks {
		prefix = mask >> shiftOffset
		if list = m[prefix]; list == nil {
			m[prefix] = &[]int{i}
		} else {
			*list = append(*list, i)
		}
	}
	lh.m5 = m

}

// IndexMasks creates a lookup table (a slice) with the p-bp prefixes,
// then you can use MaskKnownPrefixes() to masks k-mers of which the prefixes are existed.
// The lenght of prefix p can't be too big, or the lookup table (a slice) would occupy a lot of space.
func (lh *LexicHash) IndexMasks(p int) error {
	k := lh.K
	if p < 3 || p > k {
		return ErrPrefixOverflow
	}

	var prefix uint64
	var list *[]int

	m := make([]*[]int, int(math.Pow(4, float64(p))))
	shiftOffset := (k - p) << 1
	for i, mask := range lh.Masks {
		prefix = mask >> shiftOffset
		list = m[prefix]
		if list == nil {
			m[prefix] = &[]int{i}
		} else {
			*list = append(*list, i)
		}
	}
	lh.mN = m
	lh.prefix = p
	return nil
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
//     Note that k-mers of all A's or N's are skipped in k-mer generation step.
//  2. the start 0-based positions of all k-mers, with the last 1 bit as the strand
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
		rs, re = skipRegions[ri][0]-k+1, skipRegions[ri][1]
	}

	shiftOffset := (k - 3) << 1

	for {
		kmer, kmerRC, ok, _ = iter.NextKmer()
		if !ok {
			break
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
						rs, re = skipRegions[ri][0]-k+1, skipRegions[ri][1]
					}
				}

				continue
			}
		}

		if kmer == 0 || kmerRC == 0 { // all bases are A's or N's.
			continue
		}

		js = j << 1

		// ---------- positive strand ----------

		for _, i = range *m3[int(kmer>>shiftOffset)] {
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

		for _, i = range *m3[int(kmerRC>>shiftOffset)] {
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
		rs, re = skipRegions[ri][0]-k+1, skipRegions[ri][1]
	}

	for {
		kmer, kmerRC, ok, _ = iter.NextKmer()
		if !ok {
			break
		}

		if kmer == 0 || kmerRC == 0 { // all bases are A's or N's.
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
						rs, re = skipRegions[ri][0]-k+1, skipRegions[ri][1]
					}
				}
				continue
			}
		}

		js = j << 1

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

// MaskKnownPrefixes masks k-mers of which the prefixes are existed.
// So you need to run IndexMasks first.
//
// It returns
//
//  1. the list of the most similar k-mers for each mask.
//     Note that k-mers of all A's or N's are skipped in k-mer generation step.
//  2. the start 0-based positions of all k-mers, with the last 1 bit as the strand
//     flag (1 for negative strand).
//     Attention: It might be empty (len() == 0), if there's no k-mers are captured.
//
// skipRegions is optional, which is used to skip some masked regions.
// E.g., in reference indexing step, contigs of a genome can be concatenated with k-1 N's,
// where need to be ommitted.
//
// The regions should be 0-based and ascendingly sorted.
// e.g., [100, 130], [200, 230] ...
func (lh *LexicHash) MaskKnownPrefixes(s []byte, skipRegions [][2]int) (*[]uint64, *[][]int, error) {
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
	// fast locating of mask to compare with prefix indexing

	// This k-mer iterator is a simplified version of
	// https://github.com/shenwei356/bio/blob/master/sketches/iterator.go
	iter, err := iterator.NewKmerIterator(s, lh.K)
	if err != nil {
		return nil, nil, err
	}

	mN := lh.mN
	var list *[]int
	shiftOffset := (k - lh.prefix) << 1

	if checkRegion {
		ri = 0
		rs, re = skipRegions[ri][0]-k+1, skipRegions[ri][1]
	}

	for {
		kmer, kmerRC, ok, _ = iter.NextKmer()
		if !ok {
			break
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
						rs, re = skipRegions[ri][0]-k+1, skipRegions[ri][1]
					}
				}

				continue
			}
		}

		if kmer == 0 || kmerRC == 0 { // all bases are A's or N's.
			continue
		}

		js = j << 1

		// ---------- positive strand ----------

		list = mN[kmer>>shiftOffset]
		if list != nil {
			for _, i = range *list {
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
		}

		// ---------- negative strand ----------

		js |= 1 // add the strand flag to the location

		list = mN[kmerRC>>shiftOffset]
		if list != nil {
			for _, i = range *list {
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
	}

	// -----------------------------------------------------------------------------
	// some masks may not have any matches,
	// just set the k-mer to 0, then download analysis should skip these k-mers.

	for i, h := range *hashes {
		if h == math.MaxUint64 {
			(*_kmers)[i] = 0
			(*locses)[i] = (*locses)[i][:0]
		}
	}

	lh.poolHashes.Put(hashes)
	return _kmers, locses, nil
}

// MaskKmer returns the indexes of masks that possibly mask a k-mer.
// Don't forget to recycle the result via RecycleMaskKmerResult.
func (lh *LexicHash) MaskKmer(kmer uint64) *[]int {
	list := lh.poolList.Get().(*[]int)
	*list = (*list)[:0]

	var _list *[]int
	var shiftOffset int
	shiftOffset = (lh.K - lh.prefix) << 1
	if _list = lh.mN[kmer>>shiftOffset]; _list != nil {
		*list = append(*list, (*_list)...) // directly return _list is dangerous
		return list
	}

	shiftOffset = (lh.K - 5) << 1
	if _list = lh.m5[kmer>>shiftOffset]; _list != nil {
		*list = append(*list, (*_list)...) // directly return _list is dangerous
		return list
	}

	shiftOffset = (lh.K - 3) << 1
	_list = lh.m3[kmer>>shiftOffset]   // it can't be nil
	*list = append(*list, (*_list)...) // directly return _list is dangerous
	return list
}

// RecycleMaskKmerResult recycles the result of MaskKmer()
func (lh *LexicHash) RecycleMaskKmerResult(list *[]int) {
	if list != nil {
		lh.poolList.Put(list)
	}
}

// MaskLongSeqs is faster than Mask() for longer sequences by using longer 5-bp prefixes for creating the lookup table, requiring nMasks >= 1024.
func (lh *LexicHash) MaskLongSeqs(s []byte, skipRegions [][2]int) (*[]uint64, *[][]int, error) {
	if len(lh.Masks) < 1024 {
		return nil, nil, fmt.Errorf("MaskLongSeqs is not support for masks < 1024")
	}

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
	// m4 := lh.m4
	m5 := lh.m5
	var maskIdxs *[]int

	if checkRegion {
		ri = 0
		rs, re = skipRegions[ri][0]-k+1, skipRegions[ri][1]
	}

	shiftOffset5 := (k - 5) << 1
	shiftOffset3 := (k - 3) << 1

	for {
		kmer, kmerRC, ok, _ = iter.NextKmer()
		if !ok {
			break
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
						rs, re = skipRegions[ri][0]-k+1, skipRegions[ri][1]
					}
				}

				continue
			}
		}

		if kmer == 0 || kmerRC == 0 { // all bases are A's or N's.
			continue
		}

		js = j << 1

		// ---------- positive strand ----------

		maskIdxs = m5[int(kmer>>shiftOffset5)]
		if maskIdxs == nil {
			maskIdxs = m3[int(kmer>>shiftOffset3)]
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

		maskIdxs = m5[int(kmerRC>>shiftOffset5)]
		if maskIdxs == nil {
			maskIdxs = m3[int(kmerRC>>shiftOffset3)]
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
		rs, re = skipRegions[ri][0]-k+1, skipRegions[ri][1]
	}

	for {
		kmer, kmerRC, ok, _ = iter.NextKmer()
		if !ok {
			break
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
						rs, re = skipRegions[ri][0]-k+1, skipRegions[ri][1]
					}
				}
				continue
			}
		}

		if kmer == 0 || kmerRC == 0 { // all bases are A's or N's.
			continue
		}

		js = j << 1

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
