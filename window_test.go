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
	"bytes"
	"slices"
	"sync"
	"testing"

	"github.com/shenwei356/lexichash/iterator"
)

// checkWindowMask compares both outputs with the iterator and existing masking API.
func checkWindowMask(tb testing.TB, lh *LexicHash, w *WindowMasker, sequence []byte) {
	tb.Helper()
	codes, indexes, err := w.Mask(sequence)
	if err != nil {
		tb.Fatal(err)
	}
	iter, err := iterator.NewKmerIterator(sequence, lh.K)
	if err != nil {
		tb.Fatal(err)
	}
	var expectedCodes []uint64
	for {
		forward, reverse, ok, err := iter.NextKmer()
		if err != nil {
			tb.Fatal(err)
		}
		if !ok {
			break
		}
		expectedCodes = append(expectedCodes, forward, reverse)
	}
	if !slices.Equal(codes, expectedCodes) {
		tb.Fatal("window k-mer codes differ from the iterator")
	}
	selected, locations, err := lh.MaskKnownDistinctPrefixes(sequence, nil, false)
	if err != nil {
		tb.Fatal(err)
	}
	defer lh.RecycleMaskResult(selected, locations)
	expectedIndexes := make([]int, len(expectedCodes))
	for i := range expectedIndexes {
		expectedIndexes[i] = -1
	}
	for mask, positions := range *locations {
		for _, position := range positions {
			expectedIndexes[position] = mask
		}
	}
	if !slices.Equal(indexes, expectedIndexes) {
		tb.Fatal("window mask winners differ from the existing masking API")
	}
}

// TestWindowMaskMatchesDistinctPrefixes covers reused scratch, strand ties, gaps and prefix collisions.
func TestWindowMaskMatchesDistinctPrefixes(t *testing.T) {
	for _, masks := range []int{64, 20000} {
		lh, err := NewWithSeed(31, masks, 1, 0)
		if err != nil {
			t.Fatal(err)
		}
		prefix := 3
		if masks > 64 {
			prefix = 7
		}
		if err := lh.IndexMasks(prefix); err != nil {
			t.Fatal(err)
		}
		if err := lh.IndexMasksWithDistinctPrefixes(prefix + 1); err != nil {
			t.Fatal(err)
		}
		w, err := lh.NewWindowMasker()
		if err != nil {
			t.Fatal(err)
		}
		repeat := deterministicSequence(200)
		samePrefix := lh.Masks[len(lh.Masks)/2]
		// Force a collision in the direct prefix table, whose last mask wins.
		lh.Masks[0] = samePrefix ^ 1
		if err := lh.IndexMasksWithDistinctPrefixes(prefix + 1); err != nil {
			t.Fatal(err)
		}
		withGaps := deterministicSequence(2300)
		copy(withGaps[1100:], bytes.Repeat([]byte{'N'}, 100))
		sequences := [][]byte{
			deterministicSequence(10000), deterministicSequence(2300),
			bytes.Repeat(repeat, 20), reverseComplement(bytes.Repeat(repeat, 20)),
			withGaps, bytes.Repeat([]byte{'A'}, 1000), bytes.Repeat([]byte{'T'}, 1000),
			deterministicSequence(31), deterministicSequence(52), deterministicSequence(500),
		}
		for _, sequence := range sequences {
			checkWindowMask(t, lh, w, sequence)
		}
	}
}

// TestWindowMaskInvalidInput checks initialization errors and reuse after a failed window.
func TestWindowMaskInvalidInput(t *testing.T) {
	lh, err := NewWithSeed(31, 64, 1, 0)
	if err != nil {
		t.Fatal(err)
	}
	if _, err := lh.NewWindowMasker(); err == nil {
		t.Fatal("accepted unindexed masks")
	}
	if err := lh.IndexMasks(3); err != nil {
		t.Fatal(err)
	}
	if _, err := lh.NewWindowMasker(); err == nil {
		t.Fatal("accepted missing distinct-prefix index")
	}
	if err := lh.IndexMasksWithDistinctPrefixes(4); err != nil {
		t.Fatal(err)
	}
	w, err := lh.NewWindowMasker()
	if err != nil {
		t.Fatal(err)
	}
	for _, sequence := range [][]byte{nil, []byte("ACGT"), append(deterministicSequence(120), '!')} {
		if codes, indexes, err := w.Mask(sequence); err == nil || codes != nil || indexes != nil {
			t.Fatal("invalid window did not return an error and nil results")
		}
		checkWindowMask(t, lh, w, deterministicSequence(700))
	}
}

// TestWindowMaskBothStrandTies preserves repeated palindrome winners on both strands.
func TestWindowMaskBothStrandTies(t *testing.T) {
	sequence := bytes.Repeat([]byte("ACGT"), 25)
	iter, err := iterator.NewKmerIterator(sequence[:32], 32)
	if err != nil {
		t.Fatal(err)
	}
	code, reverse, _, err := iter.NextKmer()
	if err != nil || code != reverse {
		t.Fatal("invalid palindrome fixture")
	}
	iter.NextKmer() // finish and recycle the one-k-mer iterator
	masks := make([]uint64, 64)
	for i := range masks {
		masks[i] = code
	}
	lh, err := NewWithMasks(32, masks)
	if err != nil {
		t.Fatal(err)
	}
	if err := lh.IndexMasks(3); err != nil {
		t.Fatal(err)
	}
	if err := lh.IndexMasksWithDistinctPrefixes(4); err != nil {
		t.Fatal(err)
	}
	w, err := lh.NewWindowMasker()
	if err != nil {
		t.Fatal(err)
	}
	checkWindowMask(t, lh, w, sequence)
	_, indexes, err := w.Mask(sequence)
	if err != nil {
		t.Fatal(err)
	}
	for pos := 0; pos <= len(sequence)-32; pos += 4 {
		if indexes[pos*2] != len(masks)-1 || indexes[pos*2+1] != len(masks)-1 {
			t.Fatalf("lost a strand tie at position %d", pos)
		}
	}
}

// TestWindowMaskSoftMasking matches lowercase and IUPAC handling in the existing iterator.
func TestWindowMaskSoftMasking(t *testing.T) {
	previous := iterator.SupportSoftMasking
	iterator.SupportSoftMasking = true
	defer func() { iterator.SupportSoftMasking = previous }()
	lh := newLexicMapLexicHash(t)
	w, err := lh.NewWindowMasker()
	if err != nil {
		t.Fatal(err)
	}
	sequence := deterministicSequence(2300)
	copy(sequence[100:200], bytes.ToLower(sequence[100:200]))
	copy(sequence[900:], []byte("NRYSWKMBDHVnryswkmbdhv"))
	checkWindowMask(t, lh, w, sequence)
	checkWindowMask(t, lh, w, bytes.ToLower(sequence))
}

// TestWindowMaskConcurrentScratch checks independent scratch against shared immutable indexes.
func TestWindowMaskConcurrentScratch(t *testing.T) {
	lh := newLexicMapLexicHash(t)
	var wg sync.WaitGroup
	for i := 0; i < 8; i++ {
		wg.Add(1)
		go func() {
			defer wg.Done()
			w, err := lh.NewWindowMasker()
			if err != nil {
				t.Error(err)
				return
			}
			for n := 2300; n < 2500; n += 50 {
				checkWindowMask(t, lh, w, deterministicSequence(n))
			}
		}()
	}
	wg.Wait()
}

// BenchmarkWindowMask compares warmed reusable scratch with the complete legacy window workflow.
func BenchmarkWindowMask(b *testing.B) {
	lh := newLexicMapLexicHash(b)
	sequence := deterministicSequence(2300)
	b.Run("legacy", func(b *testing.B) {
		var codes []uint64
		indexes := make([]int, (len(sequence)-lh.K+1)*2)
		b.ReportAllocs()
		for n := 0; n < b.N; n++ {
			iter, err := iterator.NewKmerIterator(sequence, lh.K)
			if err != nil {
				b.Fatal(err)
			}
			codes = codes[:0]
			for {
				forward, reverse, ok, err := iter.NextKmer()
				if err != nil {
					b.Fatal(err)
				}
				if !ok {
					break
				}
				codes = append(codes, forward, reverse)
			}
			kmers, positions, err := lh.MaskKnownDistinctPrefixes(sequence, nil, false)
			if err != nil {
				b.Fatal(err)
			}
			for i := range indexes {
				indexes[i] = -1
			}
			for mask, locations := range *positions {
				for _, position := range locations {
					indexes[position] = mask
				}
			}
			lh.RecycleMaskResult(kmers, positions)
		}
	})
	b.Run("window", func(b *testing.B) {
		w, err := lh.NewWindowMasker()
		if err != nil {
			b.Fatal(err)
		}
		if _, _, err := w.Mask(sequence); err != nil {
			b.Fatal(err)
		}
		b.ReportAllocs()
		b.ResetTimer()
		for n := 0; n < b.N; n++ {
			if _, _, err := w.Mask(sequence); err != nil {
				b.Fatal(err)
			}
		}
	})
}
