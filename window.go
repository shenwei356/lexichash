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

// WindowMasker finds distinct-prefix mask winners in short sequence windows,
// such as the overlapping windows used to fill sketching deserts in LexicMap.
// Use it when the caller needs k-mer codes and their winning mask at each
// position, rather than full per-mask k-mer and location lists.
//
// Mask decodes each window once, records candidates directly by position, and
// resets only previously visited masks. This avoids initializing and traversing
// full per-mask result lists for each short window and converting those lists
// back to positions. It is most useful when a window visits few of the masks.
// It does not check shorter prefixes or accept skip regions; callers must
// exclude unwanted positions when selecting seeds from the returned candidates.
//
// Scratch is reused between calls. Each concurrent caller needs its own
// WindowMasker, while the shared LexicHash masks and indexes must remain unchanged.
type WindowMasker struct {
	lh          *LexicHash // shared masks and prefix indexes, immutable during use
	hashes      []uint64   // minimum XOR hash per mask, shared by both strands
	touched     []int      // masks to reset before processing the next window
	kmers       []uint64   // forward/reverse k-mer codes interleaved by position
	maskIndexes []int      // winning mask per interleaved position, or -1
}

// NewWindowMasker creates scratch for MaskKnownDistinctPrefixes with no shorter-
// prefix fallback. Call IndexMasks and IndexMasksWithDistinctPrefixes first.
func (lh *LexicHash) NewWindowMasker() (*WindowMasker, error) {
	if lh.mNOffsets == nil {
		return nil, fmt.Errorf("IndexMasks is not called first")
	}
	if lh.mU == nil {
		return nil, fmt.Errorf("IndexMasksWithDistinctPrefixes is not called first")
	}
	return &WindowMasker{lh: lh, hashes: slices.Clone(lh.defaultHashes)}, nil
}

// Mask returns k-mer codes and their winning masks, with position j and strand
// r stored at 2*j+r (r=0 forward, r=1 reverse). A mask index of -1 means no
// match. Results match MaskKnownDistinctPrefixes(s, nil, false), including all
// ties on both strands. Returned slices remain valid only until the next call.
func (w *WindowMasker) Mask(s []byte) ([]uint64, []int, error) {
	for _, i := range w.touched {
		w.hashes[i] = math.MaxUint64
	}
	w.touched = w.touched[:0]

	iter, err := iterator.NewKmerIterator(s, w.lh.K)
	if err != nil {
		return nil, nil, err
	}
	// Each valid k-mer start contributes one slot for each strand.
	n := (len(s) - w.lh.K + 1) * 2
	w.kmers = slices.Grow(w.kmers[:0], n)[:n]
	w.maskIndexes = slices.Grow(w.maskIndexes[:0], n)[:n]
	masks, distinct, hashes := w.lh.Masks, w.lh.mU, w.hashes
	shift := w.lh.shiftOffsetU // suffix bits to discard for the distinct-prefix lookup
	for {
		forward, reverse, ok, err := iter.NextKmer()
		if err != nil {
			return nil, nil, err
		}
		if !ok {
			break
		}
		pos := iter.Index() << 1
		w.kmers[pos], w.kmers[pos+1] = forward, reverse
		w.maskIndexes[pos], w.maskIndexes[pos+1] = -1, -1
		if forward == 0 || reverse == 0 {
			continue
		}
		// Handle each strand directly without an intermediate two-element array.
		if index := distinct[forward>>shift]; index != 0 {
			i := int(index - 1)
			w.maskIndexes[pos] = i
			hash := forward ^ masks[i]
			if hash < hashes[i] {
				// Matching prefixes make MaxUint64 an impossible XOR result,
				// so this sentinel identifies the first visit to the mask.
				if hashes[i] == math.MaxUint64 {
					w.touched = append(w.touched, i)
				}
				hashes[i] = hash
			}
		}
		if index := distinct[reverse>>shift]; index != 0 {
			i := int(index - 1)
			w.maskIndexes[pos+1] = i
			hash := reverse ^ masks[i]
			if hash < hashes[i] {
				if hashes[i] == math.MaxUint64 {
					w.touched = append(w.touched, i)
				}
				hashes[i] = hash
			}
		}
	}
	// Filter candidates against each mask's final minimum, preserving every tie.
	for pos, i := range w.maskIndexes {
		if i >= 0 && w.kmers[pos]^masks[i] != hashes[i] {
			w.maskIndexes[pos] = -1
		}
	}
	return w.kmers, w.maskIndexes, nil
}
