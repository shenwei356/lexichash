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
	"fmt"
	"slices"
	"sync"
	"testing"

	"github.com/shenwei356/lexichash/iterator"
)

// checkSketch compares compact results with the selected values of the old API.
func checkSketch(t testing.TB, lh *LexicHash, sketch *Sketcher, sequence []byte, regions []int, fallback bool) {
	t.Helper()
	got, err := sketch.Mask(sequence, regions, fallback)
	if err != nil {
		t.Fatal(err)
	}
	want, locs, err := lh.MaskKnownDistinctPrefixes(sequence, regions, fallback)
	if err != nil {
		t.Fatal(err)
	}
	defer lh.RecycleMaskResult(want, locs)
	if len(got) != len(sketch.MaskIndexes()) {
		t.Fatal("incorrect compact result length")
	}
	for slot, mask := range sketch.MaskIndexes() {
		if got[slot] != (*want)[mask] {
			t.Fatalf("mask %d: got %d, want %d", mask, got[slot], (*want)[mask])
		}
	}
}

// TestSketchMatchesDistinctPrefixes covers selections, fallback, gaps and reuse.
func TestSketchMatchesDistinctPrefixes(t *testing.T) {
	lh := newLexicMapLexicHash(t)
	selected := make([]int, 1024)
	for i := range selected {
		selected[i] = i * 19
	}
	withGaps := deterministicSequence(2300)
	copy(withGaps[1000:], bytes.Repeat([]byte{'N'}, 100))
	for _, selection := range [][]int{nil, {}, {17, 0, 19999, 3}, selected} {
		sketch, err := lh.NewSketcher(selection)
		if err != nil {
			t.Fatal(err)
		}
		for _, fallback := range []bool{false, true} {
			for _, sequence := range [][]byte{
				deterministicSequence(65536), deterministicSequence(31),
				deterministicSequence(2300), withGaps,
				bytes.Repeat([]byte{'A'}, 2300), bytes.Repeat([]byte{'T'}, 2300),
				bytes.Repeat([]byte("ACGT"), 600), deterministicSequence(52),
			} {
				checkSketch(t, lh, sketch, sequence, nil, fallback)
				if len(sequence) > 150 {
					checkSketch(t, lh, sketch, sequence, []int{0, 10, 100, 120}, fallback)
				}
			}
		}
	}
}

// TestSketchGlobalDistinctPrefix keeps fallback dependent on the complete table.
func TestSketchGlobalDistinctPrefix(t *testing.T) {
	sequence := deterministicSequence(31)
	iter, err := iterator.NewKmerIterator(sequence, 31)
	if err != nil {
		t.Fatal(err)
	}
	code, _, _, err := iter.NextKmer()
	if err != nil {
		t.Fatal(err)
	}
	iter.NextKmer()
	masks := make([]uint64, 64)
	for i := range masks {
		masks[i] = code ^ (1 << 46) // same seven-base prefix, different eight-base prefix
	}
	masks[63] = code ^ 1 // the direct match is unselected
	lh, err := NewWithMasks(31, masks)
	if err != nil {
		t.Fatal(err)
	}
	if err := lh.IndexMasks(7); err != nil {
		t.Fatal(err)
	}
	if err := lh.IndexMasksWithDistinctPrefixes(8); err != nil {
		t.Fatal(err)
	}
	// NewWithMasks sorts its input; locate the intended fallback mask afterward.
	selected := slices.Index(lh.Masks, code^(1<<46))
	sketch, err := lh.NewSketcher([]int{selected})
	if err != nil {
		t.Fatal(err)
	}
	checkSketch(t, lh, sketch, sequence, nil, true)
	got, err := sketch.Mask(sequence, nil, true)
	if err != nil || got[0] != 0 {
		t.Fatalf("fell back after an unselected direct match: got=%v err=%v", got, err)
	}
}

// TestSketchInvalidInput checks selection ownership and reuse after input errors.
func TestSketchInvalidInput(t *testing.T) {
	lh, err := NewWithSeed(31, 64, 1, 0)
	if err != nil {
		t.Fatal(err)
	}
	if _, err := lh.NewSketcher(nil); err == nil {
		t.Fatal("accepted missing prefix index")
	}
	if err := lh.IndexMasks(3); err != nil {
		t.Fatal(err)
	}
	if _, err := lh.NewSketcher(nil); err == nil {
		t.Fatal("accepted missing distinct prefix index")
	}
	if err := lh.IndexMasksWithDistinctPrefixes(4); err != nil {
		t.Fatal(err)
	}
	for _, indexes := range [][]int{{-1}, {64}, {1, 1}, {0, -1, 2}} {
		if _, err := lh.NewSketcher(indexes); err == nil {
			t.Fatal("accepted invalid selection")
		}
	}
	selection := []int{9, 1, 3}
	sketch, err := lh.NewSketcher(selection)
	if err != nil {
		t.Fatal(err)
	}
	selection[0] = 0
	if !slices.Equal(sketch.MaskIndexes(), []int{9, 1, 3}) {
		t.Fatal("selection was not copied")
	}
	checkSketch(t, lh, sketch, deterministicSequence(1000), nil, true)
	for _, sequence := range [][]byte{nil, []byte("ACGT"), append(deterministicSequence(120), '!')} {
		if result, err := sketch.Mask(sequence, nil, true); err == nil || result != nil {
			t.Fatal("invalid sequence did not return nil results and an error")
		}
		checkSketch(t, lh, sketch, deterministicSequence(700), nil, true)
	}
	if result, err := sketch.Mask(deterministicSequence(120), []int{0}, true); err != ErrInvalidSkipRegions || result != nil {
		t.Fatal("accepted an incomplete skip-region pair")
	}
	checkSketch(t, lh, sketch, deterministicSequence(700), nil, true)
}

// TestSketchSoftMasking matches lowercase and IUPAC handling in the old API.
func TestSketchSoftMasking(t *testing.T) {
	previous := iterator.SupportSoftMasking
	iterator.SupportSoftMasking = true
	defer func() { iterator.SupportSoftMasking = previous }()
	lh := newLexicMapLexicHash(t)
	sketch, err := lh.NewSketcher([]int{0, 17, 19999})
	if err != nil {
		t.Fatal(err)
	}
	sequence := deterministicSequence(2300)
	copy(sequence[100:300], bytes.ToLower(sequence[100:300]))
	copy(sequence[800:], []byte("NRYSWKMBDHVnryswkmbdhv"))
	checkSketch(t, lh, sketch, sequence, []int{900, 1000}, true)
	checkSketch(t, lh, sketch, bytes.ToLower(sequence), nil, true)
}

// TestSketchConcurrentScratch checks independent scratch against shared indexes.
func TestSketchConcurrentScratch(t *testing.T) {
	lh := newLexicMapLexicHash(t)
	var wg sync.WaitGroup
	for i := 0; i < 8; i++ {
		wg.Add(1)
		go func(mask int) {
			defer wg.Done()
			sketch, err := lh.NewSketcher([]int{mask, mask + 8})
			if err != nil {
				t.Error(err)
				return
			}
			for _, length := range []int{2300, 1542, 65536} {
				checkSketch(t, lh, sketch, deterministicSequence(length), nil, true)
			}
		}(i)
	}
	wg.Wait()
}

// BenchmarkSketcher compares full results with compact k-mer-only capture.
func BenchmarkSketcher(b *testing.B) {
	lh := newLexicMapLexicHash(b)
	for _, count := range []int{1024, 20000} {
		selection := make([]int, count)
		for i := range selection {
			selection[i] = i * len(lh.Masks) / count
		}
		for _, length := range []int{2300, 65536, 1000000} {
			sequence := deterministicSequence(length)
			b.Run(fmt.Sprintf("legacy/%d/%d", count, length), func(b *testing.B) {
				kmers, locs, err := lh.MaskKnownDistinctPrefixes(sequence, nil, true)
				if err != nil {
					b.Fatal(err)
				}
				lh.RecycleMaskResult(kmers, locs)
				b.ReportAllocs()
				b.ResetTimer()
				for n := 0; n < b.N; n++ {
					kmers, locs, err := lh.MaskKnownDistinctPrefixes(sequence, nil, true)
					if err != nil {
						b.Fatal(err)
					}
					var checksum uint64
					for _, mask := range selection {
						checksum += (*kmers)[mask]
					}
					if checksum == 0 {
						b.Fatal("empty benchmark sketch")
					}
					lh.RecycleMaskResult(kmers, locs)
				}
			})
			b.Run(fmt.Sprintf("sketch/%d/%d", count, length), func(b *testing.B) {
				sketch, err := lh.NewSketcher(selection)
				if err != nil {
					b.Fatal(err)
				}
				if _, err := sketch.Mask(sequence, nil, true); err != nil {
					b.Fatal(err)
				}
				b.ReportAllocs()
				b.ResetTimer()
				for n := 0; n < b.N; n++ {
					kmers, err := sketch.Mask(sequence, nil, true)
					if err != nil {
						b.Fatal(err)
					}
					var checksum uint64
					for _, kmer := range kmers {
						checksum += kmer
					}
					if checksum == 0 {
						b.Fatal("empty benchmark sketch")
					}
				}
			})
		}
	}
}
