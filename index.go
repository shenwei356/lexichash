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

// Index creates LexicHash index for mutitple reference sequences
// and supports searching with a query sequence.
type Index struct {
	lh *LexicHash

	// each record of the k-mer value is an uint64
	//  ref idx: 26 bits
	//  pos:     36 bits
	//  strand:   2 bits
	Trees []*tree.Tree

	IDs [][]byte // IDs of the reference genomes
	i   uint32   // curent index, for inserting a new ref seq
}

// NewIndex ceates a new Index.
func NewIndex(k int, nMasks int) (*Index, error) {
	return NewIndexWithSeed(k, nMasks, 1)
}

// NewIndexWithSeed ceates a new Index with given seed.
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

// Insert adds a new reference sequence to the index
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

// Substr represents a found substring.
type Substr struct {
	kmers.KmerCode // save substring in KmerCode

	Begin int  // start position of the substring
	End   int  // end position of the substring
	RC    bool // is the substring from the negative strand
}

// Equal tells if two Substr are the same.
func (s Substr) Equal(b Substr) bool {
	return s.K == b.K && s.Code == b.Code && s.Begin == b.Begin
}

func (s Substr) String() string {
	return fmt.Sprintf("%s %d-%d rc: %v", s.KmerCode.String(), s.Begin, s.End, s.RC)
}

// SearchResult stores a search result for the given query sequence.
type SearchResult struct {
	IdIdx int         // index of the matched reference ID
	score float64     // score for sorting
	Subs  [][2]Substr // matched substring pairs (query,target)

	scoring bool // is score computed
	cleaned bool // is duplicates removed
}

func (r SearchResult) String() string {
	return fmt.Sprintf("IdIdx: %d, Subs: %v", r.IdIdx, r.Subs)
}

// Score computes the score
func (r *SearchResult) Score() float64 {
	if r.scoring {
		return r.score
	}

	if !r.cleaned {
		r.Deduplicate()
	}

	for _, v := range r.Subs {
		r.score += float64(v[0].K)
	}

	r.scoring = true
	return r.score
}

// Deduplicate removes duplicated substrings
func (r *SearchResult) Deduplicate() {
	if len(r.Subs) == 1 {
		r.cleaned = true
		return
	}

	sort.Slice(r.Subs, func(i, j int) bool {
		if r.Subs[i][1].Begin == r.Subs[j][1].Begin {
			return r.Subs[i][1].End > r.Subs[j][1].End
		}
		return r.Subs[i][1].Begin < r.Subs[j][1].Begin
	})

	var i, j int
	var p, v [2]Substr
	var flag bool
	p = r.Subs[0]
	for i = 1; i < len(r.Subs); i++ {
		v = r.Subs[i]
		if v[1].Equal(p[1]) || v[1].End <= p[1].End { // the same or nested
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

	r.cleaned = true
}

// Search queries the index with a sequence
func (idx *Index) Search(s []byte, minPrefix uint8) ([]*SearchResult, error) {
	_kmers, _locs, err := idx.lh.Mask(s)
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

	// query substring
	var _code uint64
	var _pos int
	var _begin, _end int
	var _rc bool

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

		// queyr substring
		_pos = _locs[i]
		_rc = kmer&1 > 0

		// fmt.Printf("%3d %s\n", i, kmers.Decode(kmer>>2, k))
		for _, sr = range srs {
			// fmt.Printf("    %s %d\n",
			// 	kmers.Decode(tree.KmerPrefix(sr.Kmer, sr.K, sr.LenPrefix), int(sr.LenPrefix)),
			// 	sr.LenPrefix)

			K = int(sr.K)
			_k = int(sr.LenPrefix)

			// query
			if _rc {
				_begin, _end = _pos+K-_k, _pos+K
			} else {
				_begin, _end = _pos, _pos+_k
			}
			_code = tree.KmerPrefix(kmer>>2, sr.K, sr.LenPrefix)

			// matched
			code = tree.KmerPrefix(sr.Kmer, sr.K, sr.LenPrefix)

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
						Subs: [][2]Substr{{
							{kmers.KmerCode{Code: _code, K: _k}, _begin, _end, _rc},
							{kmers.KmerCode{Code: code, K: _k}, begin, end, rc},
						}},
					}
					m[idIdx] = r
				} else {
					r.Subs = append(r.Subs, [2]Substr{
						{kmers.KmerCode{Code: _code, K: _k}, _begin, _end, _rc},
						{kmers.KmerCode{Code: code, K: _k}, begin, end, rc},
					})
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
