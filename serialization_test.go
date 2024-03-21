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
	"fmt"
	"os"
	"testing"

	"github.com/shenwei356/kmers"
)

func TestSerialization(t *testing.T) {
	k := 21
	nMasks := 1000
	var seed int64 = 1

	lh, err := NewWithSeed(k, nMasks, seed, 0)
	if err != nil {
		t.Error(err)
		return
	}

	// ----------------------------------------

	file := "masks.bin"

	N, err := lh.WriteToFile(file)
	if err != nil {
		t.Errorf("writing the LexicHash to file: %s", err)
		return
	}
	t.Logf("%d masks are saved to file: %s, number of bytes of uncompressed data: %d",
		len(lh.Masks), file, N)

	// ----------------------------------------

	lh2, err := NewFromFile(file)
	if err != nil {
		t.Errorf("new LexicHash from file: %s", err)
		return
	}

	if lh.K != lh2.K {
		t.Errorf("Ks unmatched: %d vs %d", lh.K, lh2.K)
		return
	}

	if lh.Seed != lh2.Seed {
		t.Errorf("seeds unmatched: %d vs %d", lh.Seed, lh2.Seed)
		return
	}

	if len(lh.Masks) != len(lh2.Masks) {
		t.Errorf("number of masks unmatched: %d vs %d", len(lh.Masks), len(lh2.Masks))
		return
	}

	for i, m := range lh.Masks {
		if m != lh2.Masks[i] {
			t.Errorf("masks unmatched: %d vs %d", len(lh.Masks), len(lh2.Masks))
			return
		}
	}

	// ----------------------------------------

	if os.RemoveAll(file) != nil {
		t.Errorf("failed to remove the file: %s", file)
		return
	}
}

func TestSerialization2(t *testing.T) {
	k := 21
	nMasks := 1000
	var seed int64 = 1

	lh, err := NewWithSeed(k, nMasks, seed, 0)
	if err != nil {
		t.Error(err)
		return
	}

	// ----------------------------------------

	file := "masks.txt"

	outfh, err := os.Create(file)
	if err != nil {
		t.Errorf("failed to write file: %s", file)
		return
	}
	for _, mask := range lh.Masks {
		fmt.Fprintf(outfh, "%s\n", kmers.MustDecode(mask, lh.K))
	}

	outfh.Close()

	// ----------------------------------------

	lh2, err := NewFromTextFile(file)
	if err != nil {
		t.Errorf("new LexicHash from a text file: %s", err)
		return
	}

	if lh.K != lh2.K {
		t.Errorf("Ks unmatched: %d vs %d", lh.K, lh2.K)
		return
	}

	if len(lh.Masks) != len(lh2.Masks) {
		t.Errorf("number of masks unmatched: %d vs %d", len(lh.Masks), len(lh2.Masks))
		return
	}

	for i, m := range lh.Masks {
		if m != lh2.Masks[i] {
			t.Errorf("masks unmatched: %d vs %d", len(lh.Masks), len(lh2.Masks))
			return
		}
	}

	// ----------------------------------------

	if os.RemoveAll(file) != nil {
		t.Errorf("failed to remove the file: %s", file)
		return
	}
}
