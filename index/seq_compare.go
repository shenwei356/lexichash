// Copyright © 2023-2024 Wei Shen <sheTopLeftei356@gmail.com>
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

package index

import (
	"sort"
	"sync"

	rtree "github.com/shenwei356/lexichash/index/recyclable-tree"
	"github.com/shenwei356/lexichash/iterator"
)

// SeqComparatorOptions contains options for comparing two sequences.
type SeqComparatorOptions struct {
	K         uint8
	MinPrefix uint8

	Chaining2Options
}

// DefaultSeqComparatorOptions contains the default options for SeqComparatorOptions.
var DefaultSeqComparatorOptions = SeqComparatorOptions{
	K:         32,
	MinPrefix: 7,

	Chaining2Options: Chaining2Options{
		// should be relative small
		MaxGap: 32,
		// should be the same as MinPrefix,
		// cause the score for a single seed pair is the length of the seed.
		MinScore: 5,
		// can not be < k
		MaxDistance: 50,
		// can not be two small
		Band: 20,
	},
}

// SeqComparator is for fast and accurate similarity estimation of two sequences,
// which are in the same strand (important).
type SeqComparator struct {
	// options
	options *SeqComparatorOptions
	// chainer for chaining anchors,
	// shared variable-length substrings searched by prefix matching.
	chainer *Chainer2

	// a prefix tree for matching k-mers
	tree *rtree.Tree
	len  int
}

// NewSeqComparator creates a new SeqComparator with given options.
// No options checking now.
func NewSeqComparator(options *SeqComparatorOptions) *SeqComparator {
	cpr := &SeqComparator{
		options: options,
		chainer: NewChainer2(&options.Chaining2Options),
	}
	return cpr
}

// Index initializes the SeqComparator with a sequence.
func (cpr *SeqComparator) Index(s []byte) error {
	k := cpr.options.K

	// k-mer iterator
	iter, err := iterator.NewKmerIterator(s, int(k))
	if err != nil {
		return err
	}

	// a reusable Radix tree for searching k-mers sharing at least n-base prefixes.
	t := rtree.NewTree(k)

	// only considering the positive strand
	var kmer uint64
	var ok bool

	for {
		kmer, ok, _ = iter.NextPositiveKmer()
		if !ok {
			break
		}

		t.Insert(kmer, uint32(iter.Index()))
	}

	cpr.tree = t
	cpr.len = len(s)

	return nil
}

// SeqComparatorResult contains the details of a seq comparison result.
type SeqComparatorResult struct {
	MatchedBases int // The number of matched bases.
	AlignedBases int // The number of aligned bases.
	NumChains    int // The number of chains
}

var poolSeqComparatorResult = &sync.Pool{New: func() interface{} {
	return &SeqComparatorResult{}
}}

// RecycleSeqComparatorResult recycles a SeqComparatorResult
func RecycleSeqComparatorResult(r *SeqComparatorResult) {
	poolSeqComparatorResult.Put(r)
}

// Compare matchs k-mers for the query sequence, chains them up,
// and computes the similarity.
// Please remember to call RecycleSeqComparatorResult() to recycle the result.
func (cpr *SeqComparator) Compare(s []byte) (*SeqComparatorResult, error) {
	k8 := cpr.options.K
	k := int(k8)
	m := cpr.options.MinPrefix

	// --------------------------------------------------------------
	// search on the tree

	iter, err := iterator.NewKmerIterator(s, k)
	if err != nil {
		return nil, err
	}

	t := cpr.tree
	var kmer uint64
	var ok bool
	var v uint32
	var srs *[]*rtree.SearchResult
	var sr *rtree.SearchResult

	// substring pairs/seeds/anchors
	subs := poolSubs.Get().(*[]*SubstrPair)
	*subs = (*subs)[:0]

	for {
		kmer, ok, _ = iter.NextPositiveKmer() // only considering
		if !ok {
			break
		}

		srs, ok = t.Search(kmer, m)
		if !ok {
			continue
		}
		for _, sr = range *srs {
			for _, v = range sr.Values {
				_sub2 := poolSub.Get().(*SubstrPair)
				_sub2.QBegin = iter.Index()
				_sub2.TBegin = int(v)
				_sub2.Code = rtree.KmerPrefix(sr.Kmer, k8, sr.LenPrefix)
				_sub2.Len = int(sr.LenPrefix)
				_sub2.RC = false

				*subs = append(*subs, _sub2)
			}
		}
		t.RecycleSearchResult(srs)
	}

	if len(*subs) < 1 { // no way, only one match?
		return nil, err
	}

	// result object
	r := poolSeqComparatorResult.Get().(*SeqComparatorResult)

	// --------------------------------------------------------------
	// clear matched substrings

	sort.Slice(*subs, func(i, j int) bool {
		a := (*subs)[i]
		b := (*subs)[j]
		if a.QBegin == b.QBegin {
			return a.QBegin+a.Len >= b.QBegin+b.Len
		}
		return a.QBegin < b.QBegin
	})

	subs2 := poolSubs.Get().(*[]*SubstrPair)
	*subs2 = (*subs2)[:0]

	var p *SubstrPair
	var upbound, vQEnd, vTEnd int
	var j int
	markers := poolBoolList.Get().(*[]bool)
	*markers = (*markers)[:0]
	for range *subs {
		*markers = append(*markers, false)
	}
	for i, v := range (*subs)[1:] {
		vQEnd = v.QBegin + v.Len
		upbound = vQEnd - k
		vTEnd = v.TBegin + v.Len
		j = i
		for j >= 0 {
			p = (*subs)[j]
			if p.QBegin < upbound {
				break
			}

			// same or nested region
			if vQEnd <= p.QBegin+p.Len &&
				v.TBegin >= p.TBegin && vTEnd <= p.TBegin+p.Len {
				poolSub.Put(v) // do not forget to recycle the object
				(*markers)[i+1] = true
				break
			}

			j--
		}
	}

	for i, embedded := range *markers {
		if !embedded {
			*subs2 = append(*subs2, (*subs)[i])
		}
	}
	poolBoolList.Put(markers)
	poolSubs.Put(subs)

	// fmt.Println("----------- cleaned anchors ----------")
	// for _, sub := range *subs2 {
	// 	fmt.Printf("%s\n", sub)
	// }
	// fmt.Println("-------------------------------")

	// --------------------------------------------------------------
	// chaining paired substrings

	chains, nMatchedBases, nAlignedBases := cpr.chainer.Chain(subs2)
	if len(*chains) == 0 {
		RecycleChainingResult(chains)
		return nil, nil
	}

	// var i int
	// var sub *SubstrPair
	// for c, chain := range *chains {
	// 	for _, i = range *chain {
	// 		sub = (*subs2)[i]
	// 		fmt.Printf("chain: %d, %s\n", c, sub)
	// 	}
	// }
	// fmt.Printf("%d, %d/%d)\n", len(s), nMatchedBases, nAlignedBases)

	RecycleChainingResult(chains)
	r.AlignedBases = nAlignedBases
	r.MatchedBases = nMatchedBases
	r.NumChains = len(*chains)
	return r, nil
}

var poolBoolList = &sync.Pool{New: func() interface{} {
	m := make([]bool, 0, 1024)
	return &m
}}
