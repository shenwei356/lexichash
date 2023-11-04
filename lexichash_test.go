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
	"testing"
)

func TestHash(t *testing.T) {
	k := 31
	nMasks := 1000
	lh, err := New(k, nMasks)
	if err != nil {
		t.Error(err)
		return
	}

	// t.Logf("masks:")
	// for i, mask := range lh.Masks {
	// 	t.Logf("  %d %s %b\n", i+1, Kmer2dna(mask, k), mask)
	// }

	s1 := []byte("ATGACTGCCATGGAGGAGTCACAGTCGGATATCAGCCTCGAGCTCCCTCTGAGCCAGGAGACATTTTCAGGCTTATGGAAACTACTTCCTCCAGAAGATA")
	s2 := []byte("ATGACTGCCATGGAGGAGTCACAGcCGGATATCAGCCTtGAGCTTGAGCCAGGAGgCATTTTCAGGCTTATGGAAACTACTTCCTCCAGAcGATA")
	// s2 := []byte("TATCgTCTGGAGGAAGTAGTTTCCATAAGCCTGAAAATGcCTCCTGGCTCAAGCTCaAGGCTGATATCCGgCTGTGACTCCTCCATGGCAGTCAT")

	subs, err := lh.Compare(s1, s2, 5)
	if err != nil {
		t.Error(err)
		return
	}

	t.Logf("(location in s1) vs (location in s2) len(substring) substring")
	for _, sub := range subs {
		t.Logf("(%3d,%3d, %c) vs (%3d,%3d, %c)  %3d %s\n",
			sub[2]>>2+1, sub[2]>>2+sub[1], strands[sub[2]&1],
			sub[3]>>2+1, sub[3]>>2+sub[1], strands[sub[3]&1],
			sub[1], Kmer2dna(sub[0], int(sub[1])))
	}
}
