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
	"testing"

	"github.com/shenwei356/kmers"
)

func TestLowComplexity(t *testing.T) {
	type Case struct {
		Kmer string
		LowC bool
	}
	tests := []Case{
		{"AAAAAAA", true},
		{"CCCCCCC", true},
		{"GGGGGGG", true},
		{"TTTTTTT", true},

		{"ACAACAACAACAACA", true},
		{"ACAACAACAACAACG", true},
		{"ACAACAACAACACCG", true},
		{"ACAACAACAACACCA", true},
		{"ACAACAACAACAAAA", true},

		{"ACGACTACAGCAAAA", false},

		{"ACAAGGTACTCGCCG", false},
		{"ACAAGGTACTGGCCG", false},
		{"ACAAGGTACTCGACG", false},
		{"ACAAGGTACTATCAA", false},
		{"ACAAGGTACTCGGCG", false},
		{"ACAAGGTACTATTTT", false},

		{"ACACACCAATAGCAG", true},
	}

	var code uint64
	var r bool
	for i, test := range tests {
		code, _ = kmers.Encode([]byte(test.Kmer))
		r = IsLowComplexity(code, len(test.Kmer))
		if r != test.LowC {
			t.Errorf("[%d] %s, expected: %v, result: %v", i+1, test.Kmer, test.LowC, r)
		}
	}
}
