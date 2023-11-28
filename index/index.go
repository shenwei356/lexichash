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
	"errors"
	"fmt"
	"runtime"
	"sort"
	"sync"

	"github.com/shenwei356/lexichash"
	"github.com/shenwei356/lexichash/tree"
	"github.com/twotwotwo/sorts"
)

// ErrKConcurrentInsert occurs when calling Insert during calling BatchInsert.
var ErrKConcurrentInsert = errors.New("index: concurrent insertion")

// Index creates LexicHash index for mutitple reference sequences
// and supports searching with a query sequence.
type Index struct {
	lh *lexichash.LexicHash
	k  uint8

	// each record of the k-mer value is an uint64
	//  ref idx: 26 bits
	//  pos:     36 bits (0-based position)
	//  strand:   2 bits
	Trees []*tree.Tree

	IDs [][]byte // IDs of the reference genomes
	i   uint32   // curent index, for inserting a new ref seq

	RefSeqInfos []RefSeqInfo

	batchInsert bool
}

// GenomeInfo is a struct to store some basic information of a ref seq
type RefSeqInfo struct {
	GenomeSize int   // bases of all sequences
	Len        int   // length of contatenated sequences
	NumSeqs    int   // number of sequences
	SeqSizes   []int // sizes of sequences
}

// NewIndex ceates a new Index.
// nMasks >= 1000 is recommended.
// Setting canonicalKmer to true is recommended,
// cause it would produces more results.
func NewIndex(k int, nMasks int) (*Index, error) {
	return NewIndexWithSeed(k, nMasks, 1)
}

// NewIndexWithSeed ceates a new Index with given seed.
// nMasks >= 1000 is recommended.
// Setting canonicalKmer to true is recommended,
// cause it would produces more results.
func NewIndexWithSeed(k int, nMasks int, seed int64) (*Index, error) {
	lh, err := lexichash.NewWithSeed(k, nMasks, seed)
	if err != nil {
		return nil, err
	}

	// create a tree for each mask
	trees := make([]*tree.Tree, len(lh.Masks))
	for i := range trees {
		trees[i] = tree.New(uint8(k))
	}

	idx := &Index{
		lh:    lh,
		k:     uint8(k),
		Trees: trees,
		IDs:   make([][]byte, 0, 128),
		i:     0,

		RefSeqInfos: make([]RefSeqInfo, 0, 128),
	}

	return idx, nil
}

// K returns the K value
func (idx *Index) K() int {
	return int(idx.k)
}

// Threads is the maximum concurrency number for Insert() and ParallelizedSearch().
var Threads = runtime.NumCPU()

// Insert adds a new reference sequence to the index.
// Note that this method is not concurrency-safe,
// you can use BatchInsert, which is faster.
func (idx *Index) Insert(id []byte, s []byte, seqSize int, seqSizes []int) error {
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
		for i, kmer := range *_kmers {
			for _, loc = range (*locses)[i] {
				//  ref idx: 26 bits
				//  pos:     36 bits
				//  strand:   2 bits
				refpos = uint64(idx.i)<<38 | uint64(loc)

				idx.Trees[i].Insert(kmer, refpos)
			}
		}

		idx.IDs = append(idx.IDs, id)
		idx.RefSeqInfos = append(idx.RefSeqInfos, RefSeqInfo{
			GenomeSize: seqSize,
			Len:        seqSize + (len(seqSizes)-1)*(idx.K()-1),
			NumSeqs:    len(seqSizes),
			SeqSizes:   seqSizes,
		})
		idx.i++

		return nil
	}

	if Threads <= 0 {
		Threads = runtime.NumCPU()
	}

	var wg sync.WaitGroup
	tokens := make(chan int, Threads)

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
			var refpos uint64
			for i := start; i < end; i++ {
				kmer = (*_kmers)[i]
				for _, loc = range (*locses)[i] {
					//  ref idx: 26 bits
					//  pos:     36 bits
					//  strand:   2 bits
					refpos = uint64(idx.i)<<38 | uint64(loc)
					idx.Trees[i].Insert(kmer, refpos)
				}
			}
			wg.Done()
			<-tokens
		}(start, end)
	}
	wg.Wait()

	idx.IDs = append(idx.IDs, id)
	idx.RefSeqInfos = append(idx.RefSeqInfos, RefSeqInfo{
		GenomeSize: seqSize,
		Len:        seqSize + (len(seqSizes)-1)*(idx.K()-1),
		NumSeqs:    len(seqSizes),
		SeqSizes:   seqSizes,
	})
	idx.i++

	return nil
}

// RefSeq represents a reference sequence to insert.
type RefSeq struct {
	ID  []byte
	Seq []byte

	RefSeqSize int   // genome size
	SeqSizes   []int // lengths of each sequences
}

// MaskResult represents a mask result, it's only used in BatchInsert.
type MaskResult struct {
	ID     []byte
	Kmers  *[]uint64
	Locses *[][]int

	RefSeqSize int   // genome size
	SeqSizes   []int // lengths of each sequences
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
//	close(input)
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
			var wg sync.WaitGroup
			tokens := make(chan int, Threads)
			trees := idx.Trees
			var nMasks int
			var n int
			var j, start, end int
			var refIdx uint32
			var k int = idx.lh.K

			for m := range ch {
				nMasks = len(*(m.Kmers))
				n = nMasks/Threads + 1
				refIdx = idx.i
				for j = 0; j <= Threads; j++ {
					start, end = j*n, (j+1)*n
					if end > nMasks {
						end = nMasks
					}

					wg.Add(1)
					tokens <- 1
					go func(start, end int) {
						var kmer uint64
						var loc int
						var refpos uint64
						for i := start; i < end; i++ {
							kmer = (*m.Kmers)[i]
							for _, loc = range (*m.Locses)[i] {
								//  ref idx: 26 bits
								//  pos:     36 bits
								//  strand:   2 bits
								refpos = uint64(refIdx)<<38 | uint64(loc)
								trees[i].Insert(kmer, refpos)
							}
						}
						wg.Done()
						<-tokens
					}(start, end)
				}
				wg.Wait()

				idx.IDs = append(idx.IDs, m.ID)
				idx.RefSeqInfos = append(idx.RefSeqInfos, RefSeqInfo{
					GenomeSize: m.RefSeqSize,
					Len:        m.RefSeqSize + (len(m.SeqSizes)-1)*(k-1),
					NumSeqs:    len(m.SeqSizes),
					SeqSizes:   m.SeqSizes,
				})
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
				m.RefSeqSize = ref.RefSeqSize
				m.SeqSizes = ref.SeqSizes
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

// SubstrPair represents a pair of found substrings.
type SubstrPair struct {
	QBegin int // start position of the substring (0-based)
	QEnd   int // end position of the substring (0-based)
	TEnd   int // end position of the substring (0-based)

	QRC uint8 // a flag indicating if the substring from the negative strand (1 for yes)
	QK  uint8 // K size
	TRC uint8 // a flag indicating if the substring from the negative strand (1 for yes)
	TK  uint8 // K size

	// query
	QCode uint64 // k-mer

	// target
	TCode  uint64 // k-mer
	TBegin int    // start position of the substring (0-based)

}

// Equal tells if two SubstrPair are the same.
func (s *SubstrPair) Equal(b *SubstrPair) bool {
	return s.QK == b.QK && s.QCode == b.QCode && s.QBegin == b.QBegin && s.QRC == b.QRC &&
		s.TK == b.TK && s.TCode == b.TCode && s.TBegin == b.TBegin && s.TRC == b.TRC
}

func (s SubstrPair) String() string {
	return fmt.Sprintf("%s %d-%d rc: %d vs %s %d-%d rc: %d",
		lexichash.MustDecode(s.QCode, s.QK), s.QBegin, s.QEnd, s.QRC,
		lexichash.MustDecode(s.TCode, s.TK), s.TBegin, s.TEnd, s.TRC)
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

	r.score = 0
	for _, v := range *r.Subs {
		r.score += float64(v.QK) * float64(v.QK)
	}

	r.scoring = true
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

type pairs []*SubstrPair

func (s pairs) Len() int      { return len(s) }
func (s pairs) Swap(i, j int) { s[i], s[j] = s[j], s[i] }
func (s pairs) Less(i, j int) bool {
	a := s[i]
	b := s[j]
	if a.QBegin == b.QBegin {
		if a.QEnd == b.QEnd {
			if a.QRC == b.QRC {
				return a.TBegin < b.TBegin
			}
			return a.QRC < b.QRC
		}
		return a.QEnd > b.QEnd
	}
	return a.QBegin < b.QBegin
}

// Deduplicate removes duplicated substrings
func (r *SearchResult) Deduplicate() {
	if len(*r.Subs) == 1 {
		r.cleaned = true
		return
	}

	_subs := *r.Subs
	sort.Slice(_subs, func(i, j int) bool {
		a := _subs[i]
		b := _subs[j]
		if a.QBegin == b.QBegin {
			if a.QEnd == b.QEnd {
				if a.QRC == b.QRC {
					return a.TBegin < b.TBegin
				}
				return a.QRC < b.QRC
			}
			return a.QEnd > b.QEnd
		}
		return a.QBegin < b.QBegin
	})
	// sorts.Quicksort(pairs(*r.Subs))

	subs := poolSubs.Get().(*[]*SubstrPair)
	*subs = (*subs)[:1]

	var p, v *SubstrPair
	p = (*r.Subs)[0]
	(*subs)[0] = p
	for _, v = range (*r.Subs)[1:] {
		// if v.Equal(p) || // the same
		// 	v.QEnd <= p.QEnd { // or nested region
		if v.QEnd <= p.QEnd && v.TEnd <= p.TEnd { // same or nested region
			poolSub.Put(v)
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
	_kmers, _locses, err := idx.lh.Mask(s)
	if err != nil {
		return nil, err
	}
	defer idx.lh.RecycleMaskResult(_kmers, _locses)

	var refpos uint64
	var i int
	var kmer uint64
	// k := idx.lh.K

	// m := make(map[int]*SearchResult) // IdIdex -> result
	m := poolSearchResultsMap.Get().(*map[int]*SearchResult)
	clear(*m)

	// query substring
	var _code uint64
	var _pos int
	var _begin, _end int
	var _rc uint8

	var code uint64
	var K, _k int
	var _k8 uint8
	var idIdx, pos, begin, end int
	var rc uint8
	trees := idx.Trees
	K = idx.K()
	var locs []int
	for i, kmer = range *_kmers { // captured k-mers by the maskes
		// srs, ok := trees[i].Search(kmer, uint8(k), minPrefix)
		srs, ok := trees[i].Search(kmer, minPrefix) // each on the corresponding tree
		if !ok {
			continue
		}

		locs = (*_locses)[i]

		// fmt.Printf("%3d %s\n", i, kmers.Decode(kmer, k))
		for _, sr := range *srs { // different k-mers
			// fmt.Printf("    %s %d\n",
			// 	kmers.Decode(tree.KmerPrefix(sr.Kmer, sr.K, sr.LenPrefix), int(sr.LenPrefix)),
			// 	sr.LenPrefix)

			_k = int(sr.LenPrefix)
			_k8 = sr.LenPrefix

			// multiple locations for each QUERY k-mer,
			// but most of cases, there's only one.
			for _, _pos = range locs {
				_rc = uint8(_pos & 1)
				_pos >>= 2

				// query
				if _rc > 0 {
					_begin, _end = _pos+K-_k, _pos+K
				} else {
					_begin, _end = _pos, _pos+_k
				}

				_code = tree.KmerPrefix(kmer, uint8(K), sr.LenPrefix)

				// matched
				code = tree.KmerPrefix(sr.Kmer, uint8(K), sr.LenPrefix)

				// multiple locations for each MATCHED k-mer
				// but most of cases, there's only one.
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

					_sub2 := poolSub.Get().(*SubstrPair)
					_sub2.QCode = _code
					_sub2.QK = _k8
					_sub2.QBegin = _begin
					_sub2.QEnd = _end
					_sub2.QRC = _rc

					_sub2.TCode = code
					_sub2.TK = _k8
					_sub2.TBegin = begin
					_sub2.TEnd = end
					_sub2.TRC = rc

					var r *SearchResult
					if r, ok = (*m)[idIdx]; !ok {
						subs := poolSubs.Get().(*[]*SubstrPair)
						*subs = (*subs)[:0]

						r = poolSearchResult.Get().(*SearchResult)
						r.IdIdx = idIdx
						r.Subs = subs
						r.cleaned = false
						r.scoring = false
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
		r.Deduplicate()
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

// --------------------------------------------------------------

var poolPathResult = &sync.Pool{New: func() interface{} {
	return &Path{}
}}

var poolPathResults = &sync.Pool{New: func() interface{} {
	paths := make([]*Path, 0, 16)
	return &paths
}}

// RecyclePathResult recycles the node list.
func (idx *Index) RecyclePathResult(paths *[]*Path) {
	for _, p := range *paths {
		idx.Trees[p.TreeIdx].RecyclePathResult(p.Nodes)
		poolPathResult.Put(p)
	}
	poolPathResults.Put(paths)
}

// Path represents the path of query in a tree.
type Path struct {
	TreeIdx int
	Nodes   *[]string
	Bases   uint8
}

// Paths returned the paths in all trees.
// Do not forget to call RecyclePathResult after using the results.
func (idx *Index) Paths(key uint64, k uint8, minPrefix uint8) *[]*Path {
	var bases uint8
	// paths := make([]Path, 0, 8)
	paths := poolPathResults.Get().(*[]*Path)
	*paths = (*paths)[:0]
	for i, tree := range idx.Trees {
		var nodes *[]string
		// nodes, bases = tree.Path(key, uint8(k), minPrefix)
		nodes, bases = tree.Path(key, minPrefix)
		if bases >= minPrefix {
			// path := Path{TreeIdx: i, Nodes: nodes, Bases: bases}
			path := poolPathResult.Get().(*Path)
			path.TreeIdx = i
			path.Nodes = nodes
			path.Bases = bases
			*paths = append(*paths, path)
		}
	}
	return paths
}

// Strands could be used to output strand for a reverse complement flag
var Strands = [2]byte{'+', '-'}
