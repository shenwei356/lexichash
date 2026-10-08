// Copyright © 2023-2026 Wei Shen <shenwei356@gmail.com>
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
	"fmt"
	"math"
	"slices"

	"github.com/shenwei356/lexichash/iterator"
)

// Sketcher captures only k-mer values for a fixed selection of masks, without
// allocating or maintaining location lists. It is intended for genome screening
// and other callers that use a subset of MaskKnownDistinctPrefixes's k-mers.
// Minima are reset only for previously visited slots, and result storage is
// proportional to the selection. Each concurrent caller needs its own Sketcher;
// the shared LexicHash masks and indexes must remain unchanged during use.
type Sketcher struct {
	lh          *LexicHash // shared immutable masks and prefix indexes
	maskIndexes []int      // original mask index for each compact result slot
	slots       []int32    // original mask index -> result slot, -1 absent; nil identity
	hashes      []uint64   // minimum XOR hash per selected mask, shared by both strands
	kmers       []uint64   // selected k-mer per result slot, zero for no match
	touched     []int      // result slots whose minima need resetting before reuse
}

// NewSketcher creates reusable scratch for the given original mask indexes.
// A nil selection uses all masks; an empty non-nil selection uses none.
// Index order is preserved, duplicates and invalid indexes are rejected, and
// the selection is copied. Call both prefix-indexing methods first.
func (lh *LexicHash) NewSketcher(maskIndexes []int) (*Sketcher, error) {
	if lh.mNOffsets == nil {
		return nil, fmt.Errorf("IndexMasks is not called first")
	}
	if lh.mU == nil {
		return nil, fmt.Errorf("IndexMasksWithDistinctPrefixes is not called first")
	}
	if maskIndexes == nil {
		maskIndexes = make([]int, len(lh.Masks))
		for i := range maskIndexes {
			maskIndexes[i] = i
		}
	} else {
		maskIndexes = slices.Clone(maskIndexes)
	}
	identity := len(maskIndexes) == len(lh.Masks) // all masks in their original order
	for i, mask := range maskIndexes {
		if mask != i {
			identity = false
		}
	}
	var slots []int32
	if !identity {
		slots = make([]int32, len(lh.Masks))
		for i := range slots {
			slots[i] = -1
		}
		for slot, mask := range maskIndexes {
			if mask < 0 || mask >= len(lh.Masks) || slots[mask] >= 0 {
				return nil, fmt.Errorf("invalid or repeated mask index: %d", mask)
			}
			slots[mask] = int32(slot)
		}
	}
	hashes := make([]uint64, len(maskIndexes))
	for i := range hashes {
		hashes[i] = math.MaxUint64
	}
	return &Sketcher{lh: lh, maskIndexes: maskIndexes, slots: slots, hashes: hashes,
		kmers: make([]uint64, len(maskIndexes))}, nil
}

// MaskIndexes returns the original mask indexes corresponding to result slots.
// The returned selection is owned by the Sketcher and must not be modified.
func (s *Sketcher) MaskIndexes() []int {
	return s.maskIndexes
}

// capture updates one selected mask's minimum; unselected masks are ignored.
func (s *Sketcher) capture(mask int, code uint64) {
	slot := mask
	if s.slots != nil {
		slot = int(s.slots[mask])
		if slot < 0 {
			return
		}
	}
	hash, previous := code^s.lh.Masks[mask], s.hashes[slot]
	if hash < previous {
		if previous == math.MaxUint64 {
			s.touched = append(s.touched, slot)
		}
		s.hashes[slot], s.kmers[slot] = hash, code
	}
}

// Mask returns the selected k-mers from MaskKnownDistinctPrefixes, in the order
// reported by MaskIndexes. A zero means no match. The global distinct-prefix
// lookup decides whether shorter-prefix fallback is allowed, even when its mask
// is not selected. Skip-region and strand handling match the existing API.
// Returned storage is reused on the next call; consume or copy results first.
func (s *Sketcher) Mask(sequence []byte, skipRegions []int, checkShorterPrefix bool) ([]uint64, error) {
	if len(skipRegions)&1 != 0 {
		return nil, ErrInvalidSkipRegions
	}
	for _, slot := range s.touched {
		s.hashes[slot] = math.MaxUint64
	}
	s.touched = s.touched[:0]
	clear(s.kmers)
	iter, err := iterator.NewKmerIterator(sequence, s.lh.K)
	if err != nil {
		return nil, err
	}
	lh := s.lh
	ri, start, end := 0, 0, 0 // current inclusive skip interval in k-mer coordinates
	checkRegion := len(skipRegions) > 0
	if checkRegion {
		start, end = skipRegions[0]-lh.K+1, skipRegions[1]
	}
	for {
		forward, reverse, ok, err := iter.NextKmer()
		if err != nil {
			return nil, err
		}
		if !ok {
			break
		}
		pos := iter.Index()
		if checkRegion && start <= pos && pos <= end {
			if pos == end {
				ri += 2
				if ri == len(skipRegions) {
					checkRegion = false
				} else {
					start, end = skipRegions[ri]-lh.K+1, skipRegions[ri+1]
				}
			}
			continue
		}
		if forward == 0 || reverse == 0 {
			continue
		}
		if index := lh.mU[forward>>lh.shiftOffsetU]; index != 0 {
			s.capture(int(index-1), forward)
		} else if checkShorterPrefix {
			prefix := forward >> lh.shiftOffset
			for _, mask := range lh.mNIndexes[lh.mNOffsets[prefix]:lh.mNOffsets[prefix+1]] {
				s.capture(int(mask), forward)
			}
		}
		if index := lh.mU[reverse>>lh.shiftOffsetU]; index != 0 {
			s.capture(int(index-1), reverse)
		} else if checkShorterPrefix {
			prefix := reverse >> lh.shiftOffset
			for _, mask := range lh.mNIndexes[lh.mNOffsets[prefix]:lh.mNOffsets[prefix+1]] {
				s.capture(int(mask), reverse)
			}
		}
	}
	return s.kmers, nil
}
