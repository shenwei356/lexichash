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
	"errors"
	"fmt"
	"runtime"
	"sort"
	"sync"

	"github.com/shenwei356/kmers"
	"github.com/shenwei356/lexichash"
	"github.com/shenwei356/lexichash/tree"
)

// ErrKConcurrentInsert occurs when calling Insert during calling BatchInsert.
var ErrKConcurrentInsert = errors.New("lexichash: concurrent insertion")

// Index creates LexicHash index for mutitple reference sequences
// and supports searching with a query sequence.
type Index struct {
	lh *lexichash.LexicHash

	// each record of the k-mer value is an uint64
	//  ref idx: 26 bits
	//  pos:     36 bits
	//  strand:   2 bits
	Trees []*tree.Tree

	IDs [][]byte // IDs of the reference genomes
	i   uint32   // curent index, for inserting a new ref seq

	batchInsert bool
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
	lh, err := lexichash.NewWithSeed(k, nMasks, canonicalKmer, seed)
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

// Insert adds a new reference sequence to the index.
// Note that this method is not concurrency-safe,
// you can use BatchInsert, which is faster.
func (idx *Index) Insert(id []byte, s []byte) error {
	if idx.batchInsert {
		return ErrKConcurrentInsert
	}

	_kmers, locses, err := idx.lh.Mask(s)
	if err != nil {
		return err
	}
	defer idx.lh.RecycleMaskResult(_kmers, locses)

	if Threads == 1 {
		var loc int
		var refpos uint64
		k := uint8(idx.lh.K)
		for i, kmer := range *_kmers {
			for _, loc = range (*locses)[i] {
				//  ref idx: 26 bits
				//  pos:     36 bits
				//  strand:   2 bits
				refpos = uint64(uint64(idx.i)<<38 | uint64(loc)<<2 | kmer&1)

				idx.Trees[i].Insert(kmer>>2, k, refpos)
			}
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
				for _, loc = range (*locses)[i] {
					//  ref idx: 26 bits
					//  pos:     36 bits
					//  strand:   2 bits
					refpos := uint64(uint64(idx.i)<<38 | uint64(loc)<<2 | kmer&1)
					idx.Trees[i].Insert(kmer>>2, k, refpos)
				}
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

// RefSeq represents a reference sequence to insert.
type RefSeq struct {
	ID  []byte
	Seq []byte
}

// MaskResult represents a mask result, it's only used in BatchInsert.
type MaskResult struct {
	ID     []byte
	Kmers  *[]uint64
	Locses *[][]int
}

// BatchInsert insert a reference sequence in parallel.
// It returns:
//
//	chan RefSeq, for sending sequence.
//	sync.WaitGroup, for wait all masks being computed.
//	chan int, for waiting all the insertions to be done.
//
// Example:
//
//	input, done := BatchInsert()
//	// record is a fastx.Record//
//	_seq := make([]byte, len(record.Seq.Seq))
//	copy(_seq, record.Seq.Seq)
//	input <- RefSeq{
//		ID:  []byte(string(record.ID)),
//		Seq: _seq,
//	}
//
//	<- done
func (idx *Index) BatchInsert() (chan RefSeq, chan int) {
	if idx.batchInsert {
		panic(ErrKConcurrentInsert)
	}
	idx.batchInsert = true

	input := make(chan RefSeq, Threads)
	doneAll := make(chan int)

	poolMaskResult := &sync.Pool{New: func() interface{} {
		return &MaskResult{}
	}}

	go func() {
		ch := make(chan *MaskResult, Threads)
		doneInsert := make(chan int)

		// insert to tree
		go func() {
			for m := range ch {
				var wg sync.WaitGroup
				tokens := make(chan int, Threads)

				k := uint8(idx.lh.K)
				nMasks := len(*(m.Kmers))
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
							kmer = (*m.Kmers)[i]
							for _, loc = range (*m.Locses)[i] {
								//  ref idx: 26 bits
								//  pos:     36 bits
								//  strand:   2 bits
								refpos := uint64(uint64(idx.i)<<38 | uint64(loc)<<2 | kmer&1)
								idx.Trees[i].Insert(kmer>>2, k, refpos)
							}
						}
						wg.Done()
						<-tokens
					}(start, end)
				}
				wg.Wait()

				idx.IDs = append(idx.IDs, m.ID)
				idx.i++

				idx.lh.RecycleMaskResult(m.Kmers, m.Locses)
				poolMaskResult.Put(m)
			}
			doneInsert <- 1
		}()

		// compute mask
		var wg sync.WaitGroup
		tokens := make(chan int, Threads)

		for ref := range input {
			tokens <- 1
			wg.Add(1)
			go func(ref RefSeq) {
				_kmers, locses, err := idx.lh.Mask(ref.Seq)
				if err != nil {
					panic(err)
				}

				m := poolMaskResult.Get().(*MaskResult)
				m.Kmers = _kmers
				m.Locses = locses
				m.ID = ref.ID
				ch <- m

				wg.Done()
				<-tokens
			}(ref)
		}

		wg.Wait()
		close(ch)
		<-doneInsert

		doneAll <- 1
	}()

	// compute

	return input, doneAll
}

// ---------------------------- for Search ----------------------------------

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
	tmp := make([]*[2]Substr, 0, 1024)
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

	// sory by: start, end ,strand
	// sort.Slice(*r.Subs, func(i, j int) bool {
	// if (*r.Subs)[i][0].Begin == (*r.Subs)[j][0].Begin {
	// 	if (*r.Subs)[i][0].End == (*r.Subs)[j][0].End {
	// 		return (*r.Subs)[i][0].RC < (*r.Subs)[j][0].RC
	// 	}
	// 	return (*r.Subs)[i][0].End > (*r.Subs)[j][0].End
	// }
	// return (*r.Subs)[i][0].Begin < (*r.Subs)[j][0].Begin
	//	})

	_subs := *r.Subs
	sort.Slice(_subs, func(i, j int) bool {
		a := (_subs)[i][0]
		b := (_subs)[j][0]
		if a.Begin == b.Begin {
			if a.End == b.End {
				return a.RC < b.RC
			}
			return a.End > b.End
		}
		return a.Begin < b.Begin
	})

	subs := poolSubs.Get().(*[]*[2]Substr)
	*subs = (*subs)[:1]

	var p, v *[2]Substr
	p = (*r.Subs)[0]
	(*subs)[0] = p
	for _, v = range (*r.Subs)[1:] {
		if (v[0].Equal(p[0]) && v[1].Equal(p[1])) || // the same
			v[0].End <= p[0].End { // or nested region
			pool2Subtr.Put(v)
			continue
		}
		*subs = append(*subs, v)
		p = v
	}
	poolSubs.Put(r.Subs)
	r.Subs = subs

	r.cleaned = true
}

// RecycleSearchResult recycle search results objects
func (idx *Index) RecycleSearchResult(sr *[]*SearchResult) {
	if sr == nil {
		return
	}

	for _, r := range *sr {
		for _, v := range *r.Subs {
			pool2Subtr.Put(v)
		}
		poolSubs.Put(r.Subs)
		poolSearchResult.Put(r)
	}
	poolSearchResults.Put(sr)
}

// Search queries the index with a sequence.
// After using the result, do not forget to call RecycleSearchResult().
func (idx *Index) Search(s []byte, minPrefix uint8) (*[]*SearchResult, error) {
	_kmers, _locses, err := idx.lh.Mask(s)
	if err != nil {
		return nil, err
	}
	defer idx.lh.RecycleMaskResult(_kmers, _locses)

	var refpos uint64
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
	for i, kmer = range *_kmers {
		srs, ok := idx.Trees[i].Search(kmer>>2, uint8(k), minPrefix)
		if !ok {
			continue
		}

		_rc = uint8(kmer & 1)

		// fmt.Printf("%3d %s\n", i, kmers.Decode(kmer>>2, k))
		for _, sr := range *srs {
			// fmt.Printf("    %s %d\n",
			// 	kmers.Decode(tree.KmerPrefix(sr.Kmer, sr.K, sr.LenPrefix), int(sr.LenPrefix)),
			// 	sr.LenPrefix)

			K = int(sr.K)
			_k = int(sr.LenPrefix)

			for _, _pos = range (*_locses)[i] {
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

					_sub2 := pool2Subtr.Get().(*[2]Substr)
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

					var r *SearchResult
					if r, ok = m[idIdx]; !ok {
						subs := poolSubs.Get().(*[]*[2]Substr)
						*subs = (*subs)[:0]

						r = poolSearchResult.Get().(*SearchResult)
						r.IdIdx = idIdx
						r.Subs = subs
						r.cleaned = false
						r.scoring = false
						r.score = 0

						m[idIdx] = r
					}

					*r.Subs = append(*r.Subs, _sub2)
				}
			}
		}

		idx.Trees[i].RecycleSearchResult(srs)
	}

	if len(m) == 0 {
		return nil, nil
	}

	rs := poolSearchResults.Get().(*[]*SearchResult)
	*rs = (*rs)[:0]
	for _, r := range m {
		r.Deduplicate()
		*rs = append(*rs, r)
	}

	// sort by score, id index
	// sort.Slice(*rs, func(i, j int) bool {
	// if (*rs)[i].Score() == (*rs)[j].Score() {
	// 	return (*rs)[i].IdIdx < (*rs)[j].IdIdx
	// }
	// return (*rs)[i].Score() > (*rs)[j].Score()
	// })
	sort.Slice(*rs, func(i, j int) bool {
		a := (*rs)[i]
		b := (*rs)[j]
		if a.Score() == b.Score() {
			return a.IdIdx < b.IdIdx
		}
		return a.Score() > b.Score()
	})

	return rs, nil
}

// ---------------------------- for Debug ----------------------------------

// Path represents the path of query in a tree
type Path struct {
	TreeIdx int
	Nodes   []string
	Bases   uint8
}

// Paths returned the paths in all trees
func (idx *Index) Paths(key uint64, k uint8, minPrefix uint8) []Path {
	var bases uint8
	paths := make([]Path, 0, 8)
	for i, tree := range idx.Trees {
		var nodes []string
		nodes, bases = tree.Path(key, uint8(k), minPrefix)
		if bases >= minPrefix {
			paths = append(paths, Path{TreeIdx: i, Nodes: nodes, Bases: bases})
		}
	}
	return paths
}

// Strands could be used to output strand for a reverse complement flag
var Strands = [2]byte{'+', '-'}
