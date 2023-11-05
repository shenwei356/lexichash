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
// nMasks >= 1000 is recommended.
// Setting canonicalKmer to true is recommended,
// cause it would produces more results.
func NewIndex(k int, nMasks int, canonicalKmer bool) (*Index, error) {
	return NewIndexWithSeed(k, nMasks, canonicalKmer, 1)
}

// NewIndexWithSeed ceates a new Index with given seed.
// nMasks >= 1000 is recommended.
// Setting canonicalKmer to true is recommended,
// cause it would produces more results.
func NewIndexWithSeed(k int, nMasks int, canonicalKmer bool, seed int64) (*Index, error) {
	lh, err := NewWithSeed(k, nMasks, canonicalKmer, seed)
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

// Threads is the maximum concurrency number for Insert().
var Threads = runtime.NumCPU()

// Insert adds a new reference sequence to the index
func (idx *Index) Insert(id []byte, s []byte) error {
	_kmers, locs, err := idx.lh.Mask(s)
	if err != nil {
		return err
	}
	defer idx.lh.RecycleMaskResult(_kmers, locs)

	if Threads == 1 {
		var loc int
		var refpos uint64
		k := uint8(idx.lh.K)
		for i, kmer := range *_kmers {
			loc = (*locs)[i]

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
	nMasks := len(*_kmers)
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
				kmer = (*_kmers)[i]
				loc = (*locs)[i]

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

	Begin int   // start position of the substring
	End   int   // end position of the substring
	RC    uint8 // a flag indicating if the substring from the negative strand (1 for yes)
}

// Equal tells if two Substr are the same.
func (s Substr) Equal(b Substr) bool {
	return s.K == b.K && s.Code == b.Code && s.Begin == b.Begin && s.RC == b.RC
}

func (s Substr) String() string {
	return fmt.Sprintf("%s %d-%d rc: %d", s.KmerCode.String(), s.Begin, s.End, s.RC)
}

var pool2Subtr = &sync.Pool{New: func() interface{} {
	return &[2]Substr{}
}}

var poolSubs = &sync.Pool{New: func() interface{} {
	tmp := make([]*[2]Substr, 0, 128)
	return &tmp
}}

var poolSearchResults = &sync.Pool{New: func() interface{} {
	tmp := make([]*SearchResult, 0, 128)
	return &tmp
}}

// SearchResult stores a search result for the given query sequence.
type SearchResult struct {
	IdIdx int           // index of the matched reference ID
	score float64       // score for sorting
	Subs  *[]*[2]Substr // matched substring pairs (query,target)

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

	for _, v := range *r.Subs {
		r.score += float64(v[0].K)
	}

	r.scoring = true
	return r.score
}

// Deduplicate removes duplicated substrings
func (r *SearchResult) Deduplicate() {
	if len(*r.Subs) == 1 {
		r.cleaned = true
		return
	}

	sort.Slice(*r.Subs, func(i, j int) bool {
		if (*r.Subs)[i][0].Begin == (*r.Subs)[j][0].Begin {
			return (*r.Subs)[i][0].End > (*r.Subs)[j][0].End
		}
		return (*r.Subs)[i][0].Begin < (*r.Subs)[j][0].Begin
	})

	var i, j int
	var p, v *[2]Substr
	var flag bool
	p = (*r.Subs)[0]
	for i = 1; i < len(*r.Subs); i++ {
		v = (*r.Subs)[i]
		if (v[0].Equal(p[0]) && v[1].Equal(p[1])) || // the same
			v[0].End <= p[0].End { // or nested region
			if !flag {
				j = i // mark insertion position
				flag = true
			}
			continue
		}

		if flag { // need to insert to previous position
			pool2Subtr.Put((*r.Subs)[j])
			(*r.Subs)[j] = v
			j++
		}
		p = v
	}
	if j > 0 {
		*r.Subs = (*r.Subs)[:j]
	}

	r.cleaned = true
}

// RecycleSearchResult recycle search results objects
func (idx *Index) RecycleSearchResult(sr *[]*SearchResult) {
	var v *[2]Substr
	for _, r := range *sr {
		for _, v = range *r.Subs {
			pool2Subtr.Put(v)
		}
		poolSubs.Put(r.Subs)
	}
	poolSearchResults.Put(sr)
}

// Search queries the index with a sequence.
// After using the result, do not forget to call RecycleSearchResult()
func (idx *Index) Search(s []byte, minPrefix uint8) (*[]*SearchResult, error) {
	_kmers, _locs, err := idx.lh.Mask(s)
	if err != nil {
		return nil, err
	}
	defer idx.lh.RecycleMaskResult(_kmers, _locs)

	var srs *[]*tree.SearchResult
	var sr *tree.SearchResult
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
	var _rc uint8

	var code uint64
	var K, _k int
	var idIdx, pos, begin, end int
	var rc uint8
	var r *SearchResult
	var subs *[]*[2]Substr
	var _sub2 *[2]Substr
	for i, kmer = range *_kmers {
		srs, ok = idx.Trees[i].Search(kmer>>2, uint8(k), minPrefix)
		if !ok {
			continue
		}

		// queyr substring
		_pos = (*_locs)[i]
		_rc = uint8(kmer & 1)

		// fmt.Printf("%3d %s\n", i, kmers.Decode(kmer>>2, k))
		for _, sr = range *srs {
			// fmt.Printf("    %s %d\n",
			// 	kmers.Decode(tree.KmerPrefix(sr.Kmer, sr.K, sr.LenPrefix), int(sr.LenPrefix)),
			// 	sr.LenPrefix)

			K = int(sr.K)
			_k = int(sr.LenPrefix)

			// query
			if _rc > 0 {
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
				rc = uint8(refpos & 1)

				if rc > 0 {
					begin, end = pos+K-_k, pos+K
				} else {
					begin, end = pos, pos+_k
				}

				_sub2 = pool2Subtr.Get().(*[2]Substr)
				(*_sub2)[0].Code = _code
				(*_sub2)[0].K = _k
				(*_sub2)[0].Begin = _begin
				(*_sub2)[0].End = _end
				(*_sub2)[0].RC = _rc

				(*_sub2)[1].Code = code
				(*_sub2)[1].K = _k
				(*_sub2)[1].Begin = begin
				(*_sub2)[1].End = end
				(*_sub2)[1].RC = rc

				if r, ok = m[idIdx]; !ok {
					subs = poolSubs.Get().(*[]*[2]Substr)
					*subs = (*subs)[:0]

					r = &SearchResult{
						IdIdx: idIdx,
						// Subs: [][2]Substr{{
						// 	{kmers.KmerCode{Code: _code, K: _k}, _begin, _end, _rc},
						// 	{kmers.KmerCode{Code: code, K: _k}, begin, end, rc},
						// }},
						Subs: subs,
					}
					m[idIdx] = r
				}

				*r.Subs = append(*r.Subs, _sub2)

				// fmt.Println(r)
			}
		}

		idx.Trees[i].RecycleSearchResult(srs)
	}

	rs := poolSearchResults.Get().(*[]*SearchResult)
	*rs = (*rs)[:0]
	for _, r := range m {
		r.Deduplicate()
		*rs = append(*rs, r)
	}

	sort.Slice(*rs, func(i, j int) bool {
		return (*rs)[i].Score() > (*rs)[j].Score()
	})

	return rs, nil
}

// Path represents the path of query in a tree
type Path struct {
	TreeIdx int
	Nodes   []string
	Bases   int
}

// Paths returned the paths in all trees
func (idx *Index) Paths(key uint64, k uint8, minPrefix int) []Path {
	var bases int
	paths := make([]Path, 0, 8)
	for i, tree := range idx.Trees {
		var nodes []string
		nodes, bases = tree.Path(key, uint8(k))
		if bases >= minPrefix {
			paths = append(paths, Path{TreeIdx: i, Nodes: nodes, Bases: bases})
		}
	}
	return paths
}
