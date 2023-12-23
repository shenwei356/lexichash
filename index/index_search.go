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
	// return fmt.Sprintf("%s %d-%d vs %d-%d",
	// 	kmers.MustDecode(s.Code, s.Len),
	// 	s.QBegin+1, s.QBegin+s.Len, s.TBegin+1, s.TBegin+s.Len)
	return fmt.Sprintf("%d-%d vs %d-%d",
		s.QBegin+1, s.QBegin+s.Len, s.TBegin+1, s.TBegin+s.Len)
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
	Subs  *[]*SubstrPair // matched substring pairs (query,target)

	// more about the alignment detail
	ChainingScore float64 // chaining score
	Chains        *[]*[]int
}

func (r SearchResult) String() string {
	return fmt.Sprintf("IdIdx: %d, Subs: %v", r.IdIdx, r.Subs)
}

type SearchResults []*SearchResult

func (s SearchResults) Len() int      { return len(s) }
func (s SearchResults) Swap(i, j int) { s[i], s[j] = s[j], s[i] }
func (s SearchResults) Less(i, j int) bool {
	a := s[i]
	b := s[j]
	if a.ChainingScore == b.ChainingScore {
		if len(*a.Subs) == len(*b.Subs) {
			return a.IdIdx < b.IdIdx
		}
		return len(*a.Subs) < len(*b.Subs)
	}
	return a.ChainingScore > b.ChainingScore
}

// CleanSubstrs removes duplicated substrings pairs and computes the score.
func (r *SearchResult) CleanSubstrs() {
	if len(*r.Subs) == 1 {
		return
	}

	// sort substrings/seeds in ascending order based on the starting position
	// and in descending order based on the ending position.
	_subs := *r.Subs
	sort.Slice(_subs, func(i, j int) bool {
		a := _subs[i]
		b := _subs[j]
		if a.QBegin == b.QBegin {
			return a.QBegin+a.Len >= b.QBegin+b.Len
		}
		return a.QBegin < b.QBegin
	})

	// check all seeds pairs,
	// and remove a pair with both seeds covered by the previous one.
	subs := poolSubs.Get().(*[]*SubstrPair)
	*subs = (*subs)[:1]

	var p, v *SubstrPair
	p = (*r.Subs)[0] // the first pair
	(*subs)[0] = p

	for _, v = range (*r.Subs)[1:] {
		if v.QBegin+v.Len <= p.QBegin+p.Len {
			if v.TBegin+v.Len <= p.TBegin+p.Len { // same or nested region
				poolSub.Put(v) // do not forget to recycle the object
				continue
			}
		}

		*subs = append(*subs, v)
		p = v
	}
	poolSubs.Put(r.Subs)
	r.Subs = subs
}

// RecycleSearchResults recycles a search result object
func (idx *Index) RecycleSearchResult(r *SearchResult) {
	for _, sub := range *r.Subs {
		poolSub.Put(sub)
	}
	poolSubs.Put(r.Subs)

	for _, chain := range *r.Chains {
		poolChain.Put(chain)
	}
	poolChains.Put(r.Chains)

	poolSearchResult.Put(r)
}

// RecycleSearchResults recycles search results objects
func (idx *Index) RecycleSearchResults(sr *[]*SearchResult) {
	if sr == nil {
		return
	}

	for _, r := range *sr {
		idx.RecycleSearchResult(r)
	}
	poolSearchResults.Put(sr)
}

var poolSearchResultsMap = &sync.Pool{New: func() interface{} {
	m := make(map[int]*SearchResult, 1024)
	return &m
}}

// SetChainingOption replaces the default chaining option with a new one.
func (idx *Index) SetChainingOption(co *ChainingOptions) {
	idx.chainingOptions = co
	idx.poolChainers = &sync.Pool{New: func() interface{} {
		return NewChainer(co)
	}}
}

// SetSearchingOptions sets the searching options.
// Note that it overwrites the result of SetChainingOption.
func (idx *Index) SetSearchingOptions(so *SearchOptions) {
	idx.searchOptions = so

	co := &ChainingOptions{
		MaxGap:   so.MaxGap,
		MinScore: seedWeight(float64(so.MinSinglePrefix)),
	}
	idx.chainingOptions = co

	idx.poolChainers = &sync.Pool{New: func() interface{} {
		return NewChainer(co)
	}}
}

// SearchOptions defineds options used in searching.
type SearchOptions struct {
	// basic
	MinPrefix       uint8 // minimum prefix length, e.g, 15
	MinSinglePrefix uint8 // minimum prefix length of the single seed, e.g, 20
	TopN            int   // keep the topN scores

	// chaining
	MaxGap float64
}

// DefaultSearchOptions contains default option values.
var DefaultSearchOptions = SearchOptions{
	MinPrefix:       15,
	MinSinglePrefix: 20,
	TopN:            10,

	MaxGap: 5000,
}

// Search queries the index with a sequence.
// After using the result, do not forget to call RecycleSearchResult().
func (idx *Index) Search(s []byte) (*[]*SearchResult, error) {
	// ----------------------------------------------------------------
	// mask the query sequence
	_kmers, _locses, err := idx.lh.Mask(s, nil)
	if err != nil {
		return nil, err
	}
	defer idx.lh.RecycleMaskResult(_kmers, _locses)

	// ----------------------------------------------------------------
	// matching the captured k-mers in databases

	// a map for collecting matches for each reference: IdIdex -> result
	m := poolSearchResultsMap.Get().(*map[int]*SearchResult)
	clear(*m) // requires go >= v1.21

	var refpos uint64
	var i int
	var kmer uint64

	// query substring
	var _pos int
	var _begin int
	var _rc bool

	var code uint64
	var K, _k int
	var idIdx, pos, begin int
	trees := idx.Trees
	K = idx.K()
	K8 := uint8(K)
	var locs []int
	var srs *[]*tree.SearchResult
	var sr *tree.SearchResult
	var ok bool
	minPrefix := idx.searchOptions.MinPrefix
	for i, kmer = range *_kmers { // captured k-mers by the maskes
		srs, ok = trees[i].Search(kmer, minPrefix) // search it on the corresponding tree
		if !ok {                                   // no matcheds
			continue
		}

		locs = (*_locses)[i] // locations in the query

		// multiple locations for each QUERY k-mer,
		// but most of cases, there's only one.
		for _, _pos = range locs {
			_rc = _pos&1 > 0 // if on the reverse complement sequence
			_pos >>= 2

			// query
			if _rc { // on the negative strand
				_begin = _pos + K - _k
			} else {
				_begin = _pos
			}

			// different k-mers in subjects,
			// most of cases, there are more than one
			for _, sr = range *srs {
				// matched length
				_k = int(sr.LenPrefix)

				// matched
				code = tree.KmerPrefix(sr.Kmer, K8, sr.LenPrefix)

				// multiple locations for each MATCHED k-mer
				// but most of cases, there's only one.
				for _, refpos = range sr.Values {
					idIdx = int(refpos >> 38)
					pos = int(refpos << 26 >> 28)

					if refpos&1 > 0 {
						begin = pos + K - _k
					} else {
						begin = pos
					}

					_sub2 := poolSub.Get().(*SubstrPair)
					_sub2.QBegin = _begin
					_sub2.TBegin = begin
					_sub2.Code = code
					_sub2.Len = _k

					var r *SearchResult
					if r, ok = (*m)[idIdx]; !ok {
						subs := poolSubs.Get().(*[]*SubstrPair)
						*subs = (*subs)[:0]

						r = poolSearchResult.Get().(*SearchResult)
						r.IdIdx = idIdx
						r.Subs = subs
						r.ChainingScore = 0

						(*m)[idIdx] = r
					}

					*r.Subs = append(*r.Subs, _sub2)
				}
			}
		}

		trees[i].RecycleSearchResult(srs)
	}

	if len(*m) == 0 { // no results
		poolSearchResultsMap.Put(m)
		return nil, nil
	}

	// ----------------------------------------------------------------
	// chaining matches for all subject sequences

	minChainingScore := idx.chainingOptions.MinScore

	rs := poolSearchResults.Get().(*[]*SearchResult)
	*rs = (*rs)[:0]

	chainer := idx.poolChainers.Get().(*Chainer)
	for _, r := range *m {
		r.CleanSubstrs() // remove duplicates
		r.Chains, r.ChainingScore = chainer.Chain(r)
		if r.ChainingScore < minChainingScore {
			idx.RecycleSearchResult(r) // do not forget to recycle unused objects
			continue
		}
		*rs = append(*rs, r)
	}
	// sort subjects in descending order based on the score (simple statistics)
	sorts.Quicksort(SearchResults(*rs))

	poolSearchResultsMap.Put(m)

	// only keep the top N targets
	topN := idx.searchOptions.TopN
	if topN > 0 && len(*rs) > topN {
		var r *SearchResult
		for i := topN; i < len(*rs); i++ {
			r = (*rs)[i]

			// do not forget to recycle all related objets
			for _, sub := range *r.Subs {
				poolSub.Put(sub)
			}
			poolSubs.Put(r.Subs)
			poolSearchResult.Put(r)
		}
		*rs = (*rs)[:topN]
	}

	// // ----------------------------------------------------------------
	// // alignment
	// if idx.saveTwoBit {
	// 	var sub *SubstrPair
	// 	qlen := len(s)
	// 	var qs, qe, ts, te, begin, end int
	// 	// var q *[]byte

	// 	var paths *[]*[]int
	// 	var path *[]int

	// 	chainer := idx.poolChainers.Get().(*Chainer)

	// 	rdr := <-idx.twobitReaders

	// 	for _, r := range *rs {
	// 		// chaining
	// 		paths, _ = chainer.Chain(r)
	// 		for _, path = range *paths {

	// 			// left
	// 			sub = (*r.Subs)[(*path)[0]]
	// 			qs = sub.QBegin
	// 			ts = sub.TBegin

	// 			// right
	// 			sub = (*r.Subs)[(*path)[len(*path)-1]]
	// 			qe = sub.QBegin + sub.Len
	// 			te = sub.TBegin + sub.Len

	// 			begin = ts - qs
	// 			if begin < 0 {
	// 				begin = 0
	// 			}
	// 			end = te + qlen - qe

	// 			// q, err = rdr.SubSeq(r.IdIdx, begin, end)
	// 			// if err != nil {
	// 			// 	return rs, err
	// 			// }
	// 			fmt.Printf("subject:%s:%d-%d:%d\n", idx.IDs[r.IdIdx], begin+1, end+1, *path)
	// 		}
	// 		RecycleChainingResult(paths)
	// 	}
	// 	idx.twobitReaders <- rdr

	// }

	return rs, nil
}
