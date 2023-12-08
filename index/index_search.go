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
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OFTestSerializationTestSerialization ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.

package index

import (
	"fmt"
	"sort"
	"sync"

	"github.com/shenwei356/kmers"
	"github.com/shenwei356/lexichash/tree"
	"github.com/twotwotwo/sorts"
)

// SubstrPair represents a pair of found substrings.
type SubstrPair struct {
	QBegin int    // start position of the substring (0-based) in query
	TBegin int    // start position of the substring (0-based) in reference
	Len    int    // length
	Code   uint64 // k-mer
}

func (s SubstrPair) String() string {
	return fmt.Sprintf("%s %d-%d vs %d-%d",
		kmers.MustDecode(s.Code, s.Len),
		s.QBegin, s.QBegin+s.Len-1, s.TBegin, s.TBegin+s.Len-1)
}

var poolSub = &sync.Pool{New: func() interface{} {
	return &SubstrPair{}
}}

var poolSubs = &sync.Pool{New: func() interface{} {
	tmp := make([]*SubstrPair, 0, 128)
	return &tmp
}}

var poolSearchResult = &sync.Pool{New: func() interface{} {
	return &SearchResult{}
}}

var poolSearchResults = &sync.Pool{New: func() interface{} {
	tmp := make([]*SearchResult, 0, 128)
	return &tmp
}}

// SearchResult stores a search result for the given query sequence.
type SearchResult struct {
	IdIdx int            // index of the matched reference ID
	score float64        // score for sorting
	Subs  *[]*SubstrPair // matched substring pairs (query,target)

	UniqMatches int // because some SubstrPair result from duplication of k-mers
}

func (r SearchResult) String() string {
	return fmt.Sprintf("IdIdx: %d, Subs: %v", r.IdIdx, r.Subs)
}

// Score computes the score.
// You need to call Clean() first.
func (r *SearchResult) Score() float64 {
	return r.score
}

type SearchResults []*SearchResult

func (s SearchResults) Len() int      { return len(s) }
func (s SearchResults) Swap(i, j int) { s[i], s[j] = s[j], s[i] }
func (s SearchResults) Less(i, j int) bool {
	a := s[i]
	b := s[j]
	if a.Score() == b.Score() {
		if len(*a.Subs) == len(*b.Subs) {
			return a.IdIdx < b.IdIdx
		}
		return len(*a.Subs) < len(*b.Subs)
	}
	return a.Score() > b.Score()
}

// Clean removes duplicated substrings and compute the score.
func (r *SearchResult) Clean() {
	if len(*r.Subs) == 1 {
		r.UniqMatches = 1
		p := (*r.Subs)[0]
		r.score = float64(p.Len) * float64(p.Len)
		return
	}

	_subs := *r.Subs
	sort.Slice(_subs, func(i, j int) bool {
		a := _subs[i]
		b := _subs[j]
		if a.QBegin == b.QBegin {
			return a.QBegin+a.Len >= b.QBegin+b.Len
		}
		return a.QBegin < b.QBegin
	})

	subs := poolSubs.Get().(*[]*SubstrPair)
	*subs = (*subs)[:1]

	var p, v *SubstrPair
	p = (*r.Subs)[0]
	(*subs)[0] = p
	r.UniqMatches = 1
	r.score = float64(p.Len) * float64(p.Len)
	for _, v = range (*r.Subs)[1:] {
		if v.QBegin+v.Len <= p.QBegin+p.Len {
			if v.TBegin+v.Len <= p.TBegin+p.Len { // same or nested region
				poolSub.Put(v)
				continue
			}
		} else { // not the same query
			r.UniqMatches++

			r.score += float64(v.Len) * float64(v.Len)
		}

		*subs = append(*subs, v)
		p = v
	}
	poolSubs.Put(r.Subs)
	r.Subs = subs
}

// RecycleSearchResult recycle search results objects
func (idx *Index) RecycleSearchResult(sr *[]*SearchResult) {
	if sr == nil {
		return
	}

	for _, r := range *sr {
		for _, sub := range *r.Subs {
			poolSub.Put(sub)
		}
		poolSubs.Put(r.Subs)
		poolSearchResult.Put(r)
	}
	poolSearchResults.Put(sr)
}

var poolSearchResultsMap = &sync.Pool{New: func() interface{} {
	m := make(map[int]*SearchResult, 1024)
	return &m
}}

// Search queries the index with a sequence.
// After using the result, do not forget to call RecycleSearchResult().
func (idx *Index) Search(s []byte, minPrefix uint8) (*[]*SearchResult, error) {
	_kmers, _locses, err := idx.lh.Mask(s, nil)
	if err != nil {
		return nil, err
	}
	defer idx.lh.RecycleMaskResult(_kmers, _locses)

	var refpos uint64
	var i int
	var kmer uint64

	// m := make(map[int]*SearchResult) // IdIdex -> result
	m := poolSearchResultsMap.Get().(*map[int]*SearchResult)
	clear(*m)

	// query substring
	var _pos int
	var _begin int
	var _rc uint8

	var code uint64
	var K, _k int
	var idIdx, pos, begin int
	var rc uint8
	trees := idx.Trees
	K = idx.K()
	K8 := uint8(K)
	var locs []int
	for i, kmer = range *_kmers { // captured k-mers by the maskes
		srs, ok := trees[i].Search(kmer, minPrefix) // each on the corresponding tree
		if !ok {
			continue
		}

		locs = (*_locses)[i]

		// fmt.Printf("%3d %s\n", i, kmers.Decode(kmer, k))
		for _, sr := range *srs { // different k-mers
			_k = int(sr.LenPrefix)

			// multiple locations for each QUERY k-mer,
			// but most of cases, there's only one.
			for _, _pos = range locs {
				_rc = uint8(_pos & 1)
				_pos >>= 2

				// query
				if _rc > 0 {
					_begin = _pos + K - _k
				} else {
					_begin = _pos
				}

				// matched
				code = tree.KmerPrefix(sr.Kmer, K8, sr.LenPrefix)

				// multiple locations for each MATCHED k-mer
				// but most of cases, there's only one.
				for _, refpos = range sr.Values {
					idIdx = int(refpos >> 38)
					pos = int(refpos << 26 >> 28)
					rc = uint8(refpos & 1)

					if rc > 0 {
						begin = pos + K - _k
					} else {
						begin = pos
					}

					_sub2 := poolSub.Get().(*SubstrPair)
					_sub2.QBegin = _begin
					_sub2.Code = code
					_sub2.Len = _k
					_sub2.TBegin = begin

					var r *SearchResult
					if r, ok = (*m)[idIdx]; !ok {
						subs := poolSubs.Get().(*[]*SubstrPair)
						*subs = (*subs)[:0]

						r = poolSearchResult.Get().(*SearchResult)
						r.IdIdx = idIdx
						r.Subs = subs
						r.score = 0

						(*m)[idIdx] = r
					}

					*r.Subs = append(*r.Subs, _sub2)
				}
			}
		}

		trees[i].RecycleSearchResult(srs)
	}

	if len(*m) == 0 {
		return nil, nil
	}

	rs := poolSearchResults.Get().(*[]*SearchResult)
	*rs = (*rs)[:0]
	for _, r := range *m {
		r.Clean()

		// filter some low-confidence matches
		if r.UniqMatches == 1 && (*r.Subs)[0].Len < 20 {
			continue
		}
		*rs = append(*rs, r)
	}

	poolSearchResultsMap.Put(m)

	// sort by score, id index
	// sort.Slice(*rs, func(i, j int) bool {
	// 	a := (*rs)[i]
	// 	b := (*rs)[j]
	// 	if a.Score() == b.Score() {
	// 		return a.IdIdx < b.IdIdx
	// 	}
	// 	return a.Score() > b.Score()
	// })rs)

	sorts.Quicksort(SearchResults(*rs))

	return rs, nil
}
