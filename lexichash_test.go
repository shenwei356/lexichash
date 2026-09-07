// Copyright © 2026 Wei Shen <shenwei356@gmail.com>
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
	"slices"
	"testing"

	"github.com/shenwei356/lexichash/iterator"
)

func deterministicSequence(n int) []byte {
	bases := [...]byte{'A', 'C', 'G', 'T'}
	s := make([]byte, n)
	var x uint64 = 0x9e3779b97f4a7c15
	for i := range s {
		x = x*6364136223846793005 + 1442695040888963407
		s[i] = bases[x>>62]
	}
	return s
}

func reverseComplement(s []byte) []byte {
	rc := make([]byte, len(s))
	for i, b := range s {
		switch b {
		case 'A':
			rc[len(s)-1-i] = 'T'
		case 'C':
			rc[len(s)-1-i] = 'G'
		case 'G':
			rc[len(s)-1-i] = 'C'
		case 'T':
			rc[len(s)-1-i] = 'A'
		}
	}
	return rc
}

func newLexicMapLexicHash(tb testing.TB) *LexicHash {
	tb.Helper()

	lh, err := NewWithSeed(31, 20000, 1, 0)
	if err != nil {
		tb.Fatal(err)
	}
	if err = lh.IndexMasks(7); err != nil {
		tb.Fatal(err)
	}
	if err = lh.IndexMasksWithDistinctPrefixes(8); err != nil {
		tb.Fatal(err)
	}
	return lh
}

func TestMaskKnownDistinctPrefixesReverseComplement(t *testing.T) {
	lh := newLexicMapLexicHash(t)
	s := deterministicSequence(1542)
	rc := reverseComplement(s)

	kmers1, locses1, err := lh.MaskKnownDistinctPrefixes(s, nil, true)
	if err != nil {
		t.Fatal(err)
	}
	defer lh.RecycleMaskResult(kmers1, locses1)

	kmers2, locses2, err := lh.MaskKnownDistinctPrefixes(rc, nil, true)
	if err != nil {
		t.Fatal(err)
	}
	defer lh.RecycleMaskResult(kmers2, locses2)

	if !slices.Equal(*kmers1, *kmers2) {
		for i := range *kmers1 {
			if (*kmers1)[i] != (*kmers2)[i] {
				t.Fatalf("mask %d differs under reverse complementation: %d != %d", i, (*kmers1)[i], (*kmers2)[i])
			}
		}
	}

	for i, locs1 := range *locses1 {
		mapped := make([]int, len(locs1))
		for j, loc := range locs1 {
			pos := loc >> 1
			strand := loc & 1
			mapped[j] = ((len(s)-lh.K-pos)<<1 | (strand ^ 1))
		}
		locs2 := slices.Clone((*locses2)[i])
		slices.Sort(mapped)
		slices.Sort(locs2)
		if !slices.Equal(mapped, locs2) {
			t.Fatalf("locations of mask %d differ under reverse complementation: %v != %v", i, mapped, locs2)
		}
	}
}

func TestMaskKmerCompactPrefixIndexes(t *testing.T) {
	lh := newLexicMapLexicHash(t)
	s := deterministicSequence(256)
	iter, err := iterator.NewKmerIterator(s, lh.K)
	if err != nil {
		t.Fatal(err)
	}

	for {
		kmer, _, ok, err := iter.NextKmer()
		if err != nil {
			t.Fatal(err)
		}
		if !ok {
			break
		}

		expected := make([]int, 0, 2)
		prefixU := kmer >> lh.shiftOffsetU
		for i, mask := range lh.Masks {
			if mask>>lh.shiftOffsetU == prefixU {
				expected = append(expected[:0], i)
			}
		}
		if len(expected) == 0 {
			prefix := kmer >> lh.shiftOffset
			for i, mask := range lh.Masks {
				if mask>>lh.shiftOffset == prefix {
					expected = append(expected, i)
				}
			}
		}

		got := lh.MaskKmer(kmer)
		if !slices.Equal(*got, expected) {
			t.Fatalf("unexpected candidates for k-mer %d: %v != %v", kmer, *got, expected)
		}
		lh.RecycleMaskKmerResult(got)
	}
}

func TestMaskKnownDistinctPrefixesRecycledResults(t *testing.T) {
	lh := newLexicMapLexicHash(t)
	longSequence := deterministicSequence(1 << 16)
	shortSequence := longSequence[:1542]

	kmers, locses, err := lh.MaskKnownDistinctPrefixes(longSequence, nil, true)
	if err != nil {
		t.Fatal(err)
	}
	lh.RecycleMaskResult(kmers, locses)

	gotKmers, gotLocses, err := lh.MaskKnownDistinctPrefixes(shortSequence, nil, true)
	if err != nil {
		t.Fatal(err)
	}
	defer lh.RecycleMaskResult(gotKmers, gotLocses)

	wantLH := newLexicMapLexicHash(t)
	wantKmers, wantLocses, err := wantLH.MaskKnownDistinctPrefixes(shortSequence, nil, true)
	if err != nil {
		t.Fatal(err)
	}
	defer wantLH.RecycleMaskResult(wantKmers, wantLocses)

	if !slices.Equal(*gotKmers, *wantKmers) {
		t.Fatal("recycled k-mer results contain stale values")
	}
	for i := range *gotLocses {
		if !slices.Equal((*gotLocses)[i], (*wantLocses)[i]) {
			t.Fatalf("recycled locations for mask %d contain stale values", i)
		}
	}
}

func TestMaskSkipRegions(t *testing.T) {
	lh := newLexicMapLexicHash(t)
	sequence := deterministicSequence(1542)
	skipRegions := []int{800, 830, 100, 130}
	SortSkipRegions(skipRegions)

	tests := []struct {
		name string
		mask func([]byte, []int) (*[]uint64, *[][]int, error)
	}{
		{"MaskKnownPrefixes", lh.MaskKnownPrefixes},
		{"MaskKnownDistinctPrefixes", func(s []byte, regions []int) (*[]uint64, *[][]int, error) {
			return lh.MaskKnownDistinctPrefixes(s, regions, true)
		}},
		{"MaskKnownDistinctPrefixesWithStrandBias", func(s []byte, regions []int) (*[]uint64, *[][]int, error) {
			return lh.MaskKnownDistinctPrefixesWithStrandBias(s, regions, true)
		}},
	}

	for _, test := range tests {
		t.Run(test.name, func(t *testing.T) {
			kmers, locses, err := test.mask(sequence, skipRegions)
			if err != nil {
				t.Fatal(err)
			}
			defer lh.RecycleMaskResult(kmers, locses)

			var before, between, after bool
			for _, locations := range *locses {
				for _, location := range locations {
					start := location >> 1
					end := start + lh.K - 1
					for i := 0; i < len(skipRegions); i += 2 {
						if start <= skipRegions[i+1] && end >= skipRegions[i] {
							t.Fatalf("k-mer at [%d, %d] overlaps skipped region [%d, %d]", start, end, skipRegions[i], skipRegions[i+1])
						}
					}
					switch {
					case end < skipRegions[0]:
						before = true
					case start > skipRegions[1] && end < skipRegions[2]:
						between = true
					case start > skipRegions[3]:
						after = true
					}
				}
			}
			if !before || !between || !after {
				t.Fatalf("missing results outside skipped regions: before=%t, between=%t, after=%t", before, between, after)
			}
		})
	}
}

func TestSortSkipRegions(t *testing.T) {
	tests := []struct {
		name    string
		regions []int
		want    []int
	}{
		{"nil", nil, nil},
		{"empty", []int{}, []int{}},
		{"single", []int{100, 130}, []int{100, 130}},
		{"sorted", []int{100, 130, 400, 450, 800, 830}, []int{100, 130, 400, 450, 800, 830}},
		{"reverse", []int{800, 830, 400, 450, 100, 130}, []int{100, 130, 400, 450, 800, 830}},
		{"mixed", []int{400, 450, -20, -10, 800, 830, 100, 130}, []int{-20, -10, 100, 130, 400, 450, 800, 830}},
	}

	for _, test := range tests {
		t.Run(test.name, func(t *testing.T) {
			SortSkipRegions(test.regions)
			if !slices.Equal(test.regions, test.want) {
				t.Fatalf("unexpected regions: %v; want: %v", test.regions, test.want)
			}
		})
	}
}

func TestMaskSkipRegionsInSecondRound(t *testing.T) {
	lh, err := NewWithSeed(31, 1024, 1, 0)
	if err != nil {
		t.Fatal(err)
	}
	sequence := make([]byte, 128)
	for i := range sequence {
		sequence[i] = 'C'
	}
	skipRegions := []int{40, 70}

	tests := []struct {
		name string
		mask func([]byte, []int) (*[]uint64, *[][]int, error)
	}{
		{"Mask", lh.Mask},
		{"MaskLongSeqs", lh.MaskLongSeqs},
	}

	for _, test := range tests {
		t.Run(test.name, func(t *testing.T) {
			kmers, locses, err := test.mask(sequence, skipRegions)
			if err != nil {
				t.Fatal(err)
			}
			defer lh.RecycleMaskResult(kmers, locses)

			for i, locations := range *locses {
				if len(locations) == 0 {
					t.Fatalf("mask %d has no second-round result", i)
				}
				for _, location := range locations {
					start := location >> 1
					end := start + lh.K - 1
					if start <= skipRegions[1] && end >= skipRegions[0] {
						t.Fatalf("second-round k-mer at [%d, %d] overlaps skipped region [%d, %d]", start, end, skipRegions[0], skipRegions[1])
					}
				}
			}
		})
	}
}

func TestInvalidSkipRegions(t *testing.T) {
	lh := newLexicMapLexicHash(t)
	sequence := deterministicSequence(1542)
	tests := []struct {
		name string
		mask func([]byte, []int) (*[]uint64, *[][]int, error)
	}{
		{"Mask", lh.Mask},
		{"MaskKnownPrefixes", lh.MaskKnownPrefixes},
		{"MaskKnownDistinctPrefixes", func(s []byte, regions []int) (*[]uint64, *[][]int, error) {
			return lh.MaskKnownDistinctPrefixes(s, regions, true)
		}},
		{"MaskKnownDistinctPrefixesWithStrandBias", func(s []byte, regions []int) (*[]uint64, *[][]int, error) {
			return lh.MaskKnownDistinctPrefixesWithStrandBias(s, regions, true)
		}},
		{"MaskLongSeqs", lh.MaskLongSeqs},
	}

	for _, test := range tests {
		t.Run(test.name, func(t *testing.T) {
			kmers, locses, err := test.mask(sequence, []int{100})
			if !errors.Is(err, ErrInvalidSkipRegions) {
				t.Fatalf("expected ErrInvalidSkipRegions, got %v", err)
			}
			if kmers != nil || locses != nil {
				t.Fatal("invalid skip regions returned non-nil results")
			}
		})
	}
}

var benchmarkKmer uint64

func BenchmarkMaskKnownDistinctPrefixesLexicMap(b *testing.B) {
	for _, seqLen := range []int{1542, 1 << 20} {
		b.Run(stringSize(seqLen), func(b *testing.B) {
			lh := newLexicMapLexicHash(b)
			s := deterministicSequence(seqLen)
			b.SetBytes(int64(len(s) - lh.K + 1))
			b.ReportAllocs()
			b.ResetTimer()

			for i := 0; i < b.N; i++ {
				kmers, locses, err := lh.MaskKnownDistinctPrefixes(s, nil, true)
				if err != nil {
					b.Fatal(err)
				}
				benchmarkKmer ^= (*kmers)[i%len(*kmers)]
				lh.RecycleMaskResult(kmers, locses)
			}
		})
	}
}

func stringSize(n int) string {
	if n == 1<<20 {
		return "1MiB"
	}
	return "1542bp"
}
