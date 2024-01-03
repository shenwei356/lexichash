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
	"fmt"
	"sort"
	"sync"

	rtree "github.com/shenwei356/lexichash/index/recyclable-tree"
	"github.com/shenwei356/lexichash/iterator"
)

type SeqComparatorOptions struct {
	K         uint8
	MinPrefix uint8

	Chaining2Options
}

var DefaultSeqComparatorOptions = SeqComparatorOptions{
	K:         32,
	MinPrefix: 7,

	Chaining2Options: Chaining2Options{
		// should be relative small
		MaxGap: 32,
		// should be the same as MinPrefix,
		// cause the score fore a single seed pair is the length of the seed
		MinScore: 5,
		// can not be < k
		MaxDistance: 50,
		// can not be versy small
		Band: 20,
	},
}

type SeqComparator struct {
	options *SeqComparatorOptions
	chainer *Chainer2

	tree *rtree.Tree // prefix tree for k-mers
	len  int
}

func NewSeqComparator(options *SeqComparatorOptions) *SeqComparator {
	cpr := &SeqComparator{
		options: options,
		chainer: NewChainer2(&options.Chaining2Options),
	}
	return cpr
}

func (cpr *SeqComparator) Init(s []byte) error {
	k := cpr.options.K

	iter, err := iterator.NewKmerIterator(s, int(k))
	if err != nil {
		return err
	}

	t := rtree.NewTree(k)

	var kmer uint64 // , kmerRC
	var ok bool
	var js uint32

	for {
		kmer, _, ok, _ = iter.NextKmer()
		if !ok {
			break
		}

		js = uint32(iter.Index()) << 2
		t.Insert(kmer, js)

		// js |= 1
		// t.Insert(kmerRC, js)
	}

	cpr.tree = t
	cpr.len = len(s)

	return nil
}

func (cpr *SeqComparator) Compare(s []byte) (float64, error) {
	k8 := cpr.options.K
	k := int(k8)
	m := cpr.options.MinPrefix

	// --------------------------------------------------------------
	// search on the tree

	iter, err := iterator.NewKmerIterator(s, k)
	if err != nil {
		return 0, err
	}

	t := cpr.tree
	var kmer, code uint64 // , kmerRC
	var ok bool
	var v uint32
	var srs *[]*rtree.SearchResult
	var sr *rtree.SearchResult
	var rc bool
	var j, _k, _begin, begin int

	subs := poolSubs.Get().(*[]*SubstrPair)
	*subs = (*subs)[:0]

	for {
		kmer, _, ok, _ = iter.NextKmer()
		if !ok {
			break
		}

		j = iter.Index()

		srs, ok = t.Search(kmer, m)
		if !ok {
			continue
		}
		for _, sr = range *srs {
			for _, v = range sr.Values {
				_k = int(sr.LenPrefix)
				rc = v&1 > 0
				v >>= 2
				_begin = j // for kmer
				if rc {
					begin = int(v) + k - _k
				} else {
					begin = int(v)
				}
				code = rtree.KmerPrefix(sr.Kmer, k8, sr.LenPrefix)

				_sub2 := poolSub.Get().(*SubstrPair)
				_sub2.QBegin = _begin
				_sub2.TBegin = begin
				_sub2.Code = code
				_sub2.Len = _k
				_sub2.RC = rc

				*subs = append(*subs, _sub2)
			}
		}
		t.RecycleSearchResult(srs)

		// srs, ok = t.Search(kmerRC, m)
		// if !ok {
		// 	continue
		// }
		// for _, sr = range *srs {
		// 	for _, v = range sr.Values {
		// 		_k = int(sr.LenPrefix)
		// 		rc = v&1 > 0
		// 		v >>= 2
		// 		_begin = j + k - _k // for kmerRC
		// 		if rc {
		// 			begin = int(v) + k - _k
		// 		} else {
		// 			begin = int(v)
		// 		}
		// 		code = rtree.KmerPrefix(sr.Kmer, k8, sr.LenPrefix)

		// 		_sub2 := poolSub.Get().(*SubstrPair)
		// 		_sub2.QBegin = _begin
		// 		_sub2.TBegin = begin
		// 		_sub2.Code = code
		// 		_sub2.Len = _k
		// 		_sub2.RC = rc

		// 		*subs = append(*subs, _sub2)
		// 	}
		// }
		// t.RecycleSearchResult(srs)
	}

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

	if len(*subs) == 1 { // no need to clean
		*subs2 = append(*subs2, (*subs)[0])
	} else {
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
	}

	// fmt.Println("----------- cleaned anchors ----------")
	// for _, sub := range *subs2 {
	// 	fmt.Printf("%s\n", sub)
	// }
	// fmt.Println("-------------------------------")

	// --------------------------------------------------------------
	// chaining paired substrings

	chains := cpr.chainer.Chain(subs2)
	if len(*chains) == 0 {
		RecycleChainingResult(chains)
		return 0, nil
	}

	// --------------------------------------------------------------
	// compute similarity

	var matches = 0

	var i int
	var sub *SubstrPair
	for c, chain := range *chains {
		for _, i = range *chain {
			sub = (*subs2)[i]
			fmt.Printf("chain: %d, %s\n", c, sub)
		}
	}

	var ident float64
	if len(s) < cpr.len {
		ident = float64(matches) / float64(len(s))
	} else {
		ident = float64(matches) / float64(cpr.len)
	}

	RecycleChainingResult(chains)
	return ident, nil
}

var poolBoolList = &sync.Pool{New: func() interface{} {
	m := make([]bool, 0, 1024)
	return &m
}}
