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
	"fmt"

	"github.com/shenwei356/kmers"
	tree "github.com/shenwei356/lexichash/kmer-radix-tree"
)

type Index struct {
	lh *LexicHash

	// each record is an uint64
	//  ref idx: 26 bits
	//  pos:     36 bits
	//  strand:   2 bits
	Trees []*tree.Tree

	ids [][]byte // IDs
	i   uint32   // curent index
}

func NewIndex(k int, nMasks int) (*Index, error) {
	return NewIndexWithSeed(k, nMasks, 1)
}

func NewIndexWithSeed(k int, nMasks int, seed int64) (*Index, error) {
	lh, err := NewWithSeed(k, nMasks, seed)
	if err != nil {
		return nil, err
	}

	// create a tree for each mask
	trees := make([]*tree.Tree, len(lh.Masks))
	for i := range trees {
		trees[i] = tree.New()
	}

	idx := &Index{
		lh:    lh,
		Trees: trees,
		ids:   make([][]byte, 0, 128),
		i:     0,
	}
	return idx, nil
}

func (idx *Index) Insert(id []byte, s []byte) error {
	_kmers, locs, err := idx.lh.Mask(s)
	if err != nil {
		return err
	}

	var loc int
	var refpos uint64
	k := uint8(idx.lh.K)
	for i, kmer := range _kmers {
		loc = locs[i]

		//  ref idx: 26 bits
		//  pos:     36 bits
		//  strand:   2 bits
		refpos = uint64(uint64(idx.i)<<38 | uint64(loc)<<2 | kmer&1)

		idx.Trees[i].Insert(kmer>>2, k, refpos)
	}

	idx.ids = append(idx.ids, id)
	idx.i++

	return nil
}

func (idx *Index) Search(s []byte, m uint8) error {
	_kmers, _, err := idx.lh.Mask(s)
	if err != nil {
		return err
	}

	var srs []tree.SearchResult
	var sr tree.SearchResult
	var refpos uint64
	var ok bool
	var i int
	var kmer uint64
	k := idx.lh.K

	for i, kmer = range _kmers {
		srs, ok = idx.Trees[i].Search(kmer>>2, uint8(k), m)
		if !ok {
			continue
		}

		fmt.Printf("%3d %s\n", i, kmers.Decode(kmer>>2, k))
		for _, sr = range srs {
			fmt.Printf("    %s %d\n",
				kmers.Decode(tree.KmerPrefix(sr.Kmer, sr.K, sr.LenPrefix), int(sr.LenPrefix)),
				sr.LenPrefix)
			for _, refpos = range sr.Values {
				fmt.Printf("      %s, %d, %c\n",
					idx.ids[refpos>>38], refpos<<26>>28, strands[refpos&1])
			}
		}
	}

	return nil
}
