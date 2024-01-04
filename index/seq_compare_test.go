// Copyright © 2023-2024 Wei Shen <sheTopLeftei356@gmail.com>
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

package index

import (
	"testing"

	"github.com/shenwei356/lexichash/index/align"
)

func TestSeqCompare(t *testing.T) {
	// Identities: 271/288(94%)

	// Query  8    AGGTCCTGCCCCGCGACCTGCACGCCGAATACGTAGCGGCGATCGCCTTAGTCGGTACAG  67
	//             |||||||||||||||||||||||||| ||||||||   |||||||||||||||  |||||
	// Sbjct  15   AGGTCCTGCCCCGCGACCTGCACGCC-AATACGTA-TAGCGATCGCCTTAGTC--TACAG  70

	// Query  68   CCCTGGAAAACATGGCCACCGAAGTTCGTTCCCTGCAACGGACCGAAATCCACGAAGTCG  127
	//             |||||||||||||||||||||||||||||| ||||||||||||||||||||||  |||||
	// Sbjct  71   CCCTGGAAAACATGGCCACCGAAGTTCGTT-CCTGCAACGGACCGAAATCCACTGAGTCG  129

	// Query  128  AAGAACACTTTGCTAAGGGCCAAAAGGGCTCGTCAGCCATGCCGCACAAGCGGAACCCAA  187
	//                || |||| |||||||||||||||||||||||||||||||||||||||||||||||||
	// Sbjct  130  --CAATACTTCGCTAAGGGCCAAAAGGGCTCGTCAGCCATGCCGCACAAGCGGAACCCAA  187

	// Query  188  TTGGCTCCGAAAACATCTGCGGCTGTGCCCGGGTCCTGCGGGGCAACGTGGTGACCGCCT  247
	//             ||||||||||||||||||||||||||||||||||||||||||| ||||||||||||||||
	// Sbjct  188  TTGGCTCCGAAAACATCTGCGGCTGTGCCCGGGTCCTGCGGGG-AACGTGGTGACCGCCT  246

	// Query  248  ACGAAGACGTGACCCTCTGGCACGAACGCGACATCTCCCACTCCAGTG  295
	//             ||||||||||||||||  ||||||||||||||||||||||||||||||
	// Sbjct  247  ACGAAGACGTGACCCTTCGGCACGAACGCGACATCTCCCACTCCAGTG  294

	s1 := []byte("GGTTACGTATTGCTAGGTCCTGCCCCGCGACCTGCACGCCAATACGTATAGCGATCGCCTTAGTCTACAGCCCTGGAAAACATGGCCACCGAAGTTCGTTCCTGCAACGGACCGAAATCCACTGAGTCGCAATACTTCGCTAAGGGCCAAAAGGGCTCGTCAGCCATGCCGCACAAGCGGAACCCAATTGGCTCCGAAAACATCTGCGGCTGTGCCCGGGTCCTGCGGGGAACGTGGTGACCGCCTACGAAGACGTGACCCTTCGGCACGAACGCGACATCTCCCACTCCAGTGAGCAATACGTAACTGAACGAAGAACATCCGCAAAAAAAA")
	s2 := []byte("TCCACCCAGGTCCTGCCCCGCGACCTGCACGCCGAATACGTAGCGGCGATCGCCTTAGTCGGTACAGCCCTGGAAAACATGGCCACCGAAGTTCGTTCCCTGCAACGGACCGAAATCCACGAAGTCGAAGAACACTTTGCTAAGGGCCAAAAGGGCTCGTCAGCCATGCCGCACAAGCGGAACCCAATTGGCTCCGAAAACATCTGCGGCTGTGCCCGGGTCCTGCGGGGCAACGTGGTGACCGCCTACGAAGACGTGACCCTCTGGCACGAACGCGACATCTCCCACTCCAGTGCCGAACGGATGATTCTGCCGGACTCCACGGCGCTGTTG")

	// alignment
	alg := align.NewAligner(&align.AlignOptions{
		MatchScore:     1,
		MisMatchScore:  -1,
		GapScore:       -1,
		SaveAlignments: true,
		SaveMatrix:     false,
	})
	r := alg.Global(s2, s1)

	t.Logf("matches: %d, gaps: %d, len: %d, identity: %.2f%%\n",
		r.Matches, r.Gaps, r.Len, float64(r.Matches)/float64(r.Len)*100)

	// compare

	cpr := NewSeqComparator(&DefaultSeqComparatorOptions)

	err := cpr.Index(s1)
	if err != nil {
		t.Logf("%s", err)
		return
	}

	cr, err := cpr.Compare(s2)
	if err != nil {
		t.Logf("%s", err)
		return
	}
	t.Logf("nChains: %d, matched bases: %d, aligned bases: %d\n", cr.NumChains, cr.MatchedBases, cr.AlignedBases)
	RecycleSeqComparatorResult(cr)
}

//
