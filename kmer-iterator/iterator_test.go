// Copyright © 2018-2021 Wei Shen <shenwei356@gmail.com>
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

package iterater

import (
	"testing"
)

func TestKmerIterator(t *testing.T) {
	_s := "AAGTTTGAATCATTCAACTATCTAGTTTTCAGAGAACAATGTTCTCTAAAGAATAGAAAAGAGTCATTGTGCGGTGATGATGGCGGGAAGGATCCACCTG"
	sequence := []byte(_s)
	k := 10

	canonical := true
	iter, err := NewKmerIterator(sequence, k, canonical)
	if err != nil {
		t.Errorf("fail to create aa iter rator")
	}

	var _codes [2]uint64
	var ok bool
	// var idx int
	codes := make([]uint64, 0, 1024)
	for {
		_codes, ok, err = iter.NextKmer()
		if err != nil {
			t.Error(err)
		}
		if !ok {
			break
		}

		// idx = iter.Index()
		// fmt.Printf("kmer: %d-%s, %s-%b, RC:%v\n",
		// idx, iter.s.Seq[idx:idx+k], kmers.Decode(code>>2, k), code>>2, code&1 > 0)

		codes = append(codes, _codes[0])
		if !canonical {
			codes = append(codes, _codes[1])
		}
	}

	if canonical {
		if len(codes) != len(_s)-k+1 {
			t.Errorf("k-mers number error")
		}
	} else if len(codes) != (len(_s)-k+1)*2 {
		t.Errorf("k-mers number error")
	}
}
