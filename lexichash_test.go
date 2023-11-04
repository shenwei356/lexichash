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
	"strconv"
	"testing"
)

func TestHash(t *testing.T) {
	k := 21
	nMasks := 1000
	cannonical := false
	var seed int64 = 1

	lh, err := NewWithSeed(k, nMasks, cannonical, seed)
	if err != nil {
		t.Error(err)
		return
	}

	// t.Logf("masks:")
	// var d int
	// var p uint64
	// for i, mask := range lh.Masks {
	// 	d = bits.OnesCount64(mask ^ p)
	// 	t.Logf("  %3d %s %2d %42b\n", i+1, Kmer2dna(mask, k), d, mask)
	// 	p = mask
	// }

	s1 := []byte("AGAAGGACGTGGACGTGGATGCCGATAAGAAGGAGCCGTAAGGTACCGGGCGTGGGGAGGGCAGGGGCAGGGACGGGGATCAGGGGCAGCTGATCCCCGT")
	s2 := []byte("AGAAGGACGTGGACGTGGATcCCGATAAGAAGGAcGCCGTAAGGTACCaGGCGTGGGGAGGGCAGGGGaAGGGACGGGGATCAGGGGCAGaTGATCCCCGT")

	minLen := 13
	subs, err := lh.Compare(s2, s1, minLen)
	if err != nil {
		t.Error(err)
		return
	}

	t.Logf("(location in s1) vs (location in s2) len(substring) substring")
	for _, sub := range subs {
		t.Logf("(%3d,%3d, %c) vs (%3d,%3d, %c)  %3d %s\n",
			sub[2]>>2+1, sub[2]>>2+sub[1], Strands[sub[2]&1],
			sub[3]>>2+1, sub[3]>>2+sub[1], Strands[sub[3]&1],
			sub[1], Kmer2dna(sub[0], int(sub[1])))
	}

	// use the same sequence to build the index

	idx, err := NewIndexWithSeed(k, nMasks, cannonical, seed)
	if err != nil {
		t.Error(err)
		return
	}
	idx.Insert([]byte("s1"), s1)

	sr, err := idx.Search(s2, uint8(minLen))
	if err != nil {
		t.Log(err)
		return
	}
	t.Log()
	t.Logf("query: %s", "s2")
	for i, r := range sr {
		t.Logf("%4s %s\n", "#"+strconv.Itoa(i+1), idx.IDs[r.IdIdx])
		for _, v := range r.Subs[0:] {
			t.Logf("     (%3d,%3d, %c) vs (%3d,%3d, %c) %3d %s\n",
				v[0].Begin+1, v[0].End, Strands[v[0].RC],
				v[1].Begin+1, v[1].End, Strands[v[1].RC],
				v[0].K, v[0].KmerCode)
		}
	}
}
