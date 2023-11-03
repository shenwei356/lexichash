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
	"runtime"
	"sort"
	"sync"

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

	IDs [][]byte // IDs
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
		IDs:   make([][]byte, 0, 128),
		i:     0,
	}
	return idx, nil
}

func (idx *Index) Insert(id []byte, s []byte) error {
	_kmers, locs, err := idx.lh.Mask(s)
	if err != nil {
		return err
	}

	if Threads == 1 {
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

		idx.IDs = append(idx.IDs, id)
		idx.i++

		return nil
	}

	if Threads <= 0 {
		Threads = runtime.NumCPU()
	}

	var wg sync.WaitGroup
	tokens := make(chan int, Threads)

	k := uint8(idx.lh.K)
	nMasks := len(_kmers)
	n := nMasks/Threads + 1
	var start, end int
	for j := 0; j <= Threads; j++ {
		start, end = j*n, (j+1)*n
		if end > nMasks {
			end = nMasks
		}

		wg.Add(1)
		tokens <- 1
		go func(start, end int) {
			var kmer uint64
			var loc int
			for i := start; i < end; i++ {
				kmer = _kmers[i]
				loc = locs[i]

				//  ref idx: 26 bits
				//  pos:     36 bits
				//  strand:   2 bits
				refpos := uint64(uint64(idx.i)<<38 | uint64(loc)<<2 | kmer&1)
				idx.Trees[i].Insert(kmer>>2, k, refpos)
			}
			wg.Done()
			<-tokens
		}(start, end)
	}
	wg.Wait()

	idx.IDs = append(idx.IDs, id)
	idx.i++

	return nil
}

type Substr struct {
	kmers.KmerCode // save substring in KmerCode

	Begin int // start position of the substring
	End   int
	RC    bool // is the substring from the negative strand
}

func (s Substr) Equal(b Substr) bool {
	return s.K == b.K && s.Code == b.Code && s.Begin == b.Begin
}

func (s Substr) String() string {
	return fmt.Sprintf("%s, pos: %d-%d, rc: %v", s.KmerCode.String(), s.Begin, s.End, s.RC)
}

type SearchResult struct {
	IdIdx int
	score float64
	Subs  []Substr

	scoring bool
	cleaned bool
}

func (r SearchResult) String() string {
	return fmt.Sprintf("IdIdx: %d, Subs: %v", r.IdIdx, r.Subs)
}

func (r *SearchResult) Score() float64 {
	if r.scoring {
		return r.score
	}

	if !r.cleaned {
		r.Deduplicate()
	}

	// compute score

	r.scoring = true
	return r.score
}

func (r *SearchResult) Deduplicate() {
	if len(r.Subs) == 1 {
		r.cleaned = true
		return
	}

	sort.Slice(r.Subs, func(i, j int) bool {
		if r.Subs[i].Begin == r.Subs[j].Begin {
			return r.Subs[i].End > r.Subs[j].End
		}
		return r.Subs[i].Begin < r.Subs[j].Begin
	})

	var i, j int
	var p, v Substr
	var flag bool
	p = r.Subs[0]
	for i = 1; i < len(r.Subs); i++ {
		v = r.Subs[i]
		if v.Equal(p) || v.End <= p.End { // the same or nested
			if !flag {
				j = i // mark insertion position
				flag = true
			}
			continue
		}

		if flag { // need to insert to previous position
			r.Subs[j] = v
			j++
		}
		p = v
	}
	if j > 0 {
		r.Subs = r.Subs[:j]
	}

	for _, v := range r.Subs {
		r.score += float64(v.K)
	}

	r.cleaned = true
}

func (idx *Index) Search(s []byte, minPrefix uint8) ([]*SearchResult, error) {
	_kmers, _, err := idx.lh.Mask(s)
	if err != nil {
		return nil, err
	}

	var srs []tree.SearchResult
	var sr tree.SearchResult
	var refpos uint64
	var ok bool
	var i int
	var kmer uint64
	k := idx.lh.K

	m := make(map[int]*SearchResult) // IdIdex -> result

	var code uint64
	var K, _k int
	var idIdx, pos, begin, end int
	var rc bool
	var r *SearchResult
	for i, kmer = range _kmers {
		srs, ok = idx.Trees[i].Search(kmer>>2, uint8(k), minPrefix)
		if !ok {
			continue
		}

		// fmt.Printf("%3d %s\n", i, kmers.Decode(kmer>>2, k))
		for _, sr = range srs {
			// fmt.Printf("    %s %d\n",
			// 	kmers.Decode(tree.KmerPrefix(sr.Kmer, sr.K, sr.LenPrefix), int(sr.LenPrefix)),
			// 	sr.LenPrefix)

			code = tree.KmerPrefix(sr.Kmer, sr.K, sr.LenPrefix)
			K = int(sr.K)
			_k = int(sr.LenPrefix)

			for _, refpos = range sr.Values {
				// fmt.Printf("      %s, %d, %c\n",
				// 	idx.ids[refpos>>38], refpos<<26>>28, strands[refpos&1])

				idIdx = int(refpos >> 38)
				pos = int(refpos << 26 >> 28)
				rc = refpos&1 > 0

				if rc {
					begin, end = pos+K-_k, pos+K
				} else {
					begin, end = pos, pos+_k
				}

				if r, ok = m[idIdx]; !ok {
					r = &SearchResult{
						IdIdx: idIdx,
						Subs:  []Substr{{kmers.KmerCode{Code: code, K: _k}, begin, end, rc}},
					}
					m[idIdx] = r
				} else {
					r.Subs = append(r.Subs, Substr{kmers.KmerCode{Code: code, K: _k}, begin, end, rc})
				}
				// fmt.Println(r)
			}
		}
	}

	rs := make([]*SearchResult, 0, len(m))
	for _, r := range m {
		r.Deduplicate()
		rs = append(rs, r)
	}
	sort.Slice(rs, func(i, j int) bool {
		return rs[i].Score() > rs[j].Score()
	})

	return rs, nil
}
