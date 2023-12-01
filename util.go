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
	"github.com/shenwei356/kmers"
)

func Kmer2dna(code uint64, k int) []byte {
	return kmers.Decode(code, k)
}

var bit2base = [4]byte{'A', 'C', 'G', 'T'}

// MustDecoder returns a Decode function, which reuses the byte slice
func MustDecoder() func(code uint64, k uint8) []byte {
	buf := make([]byte, 32)

	return func(code uint64, k uint8) []byte {
		kmer := buf[:k]
		var i uint8
		for i = 0; i < k; i++ {
			kmer[k-1-i] = bit2base[code&3]

			// it's slower than bit2base[code&3], according to the test.
			// kmer[k-1-i] = byte(uint64(1413956417>>((code&3)<<3)) & 255)

			code >>= 2
		}
		return kmer
	}
}

// MustDecode return k-mer string
func MustDecode(code uint64, k uint8) []byte {
	kmer := make([]byte, k)
	var i uint8
	for i = 0; i < k; i++ {
		kmer[k-1-i] = bit2base[code&3]

		// it's slower than bit2base[code&3], according to the test.
		// kmer[k-1-i] = byte(uint64(1413956417>>((code&3)<<3)) & 255)

		code >>= 2
	}
	return kmer
}

// Strands could be used to output strand for a reverse complement flag
var Strands = [2]byte{'+', '-'}
