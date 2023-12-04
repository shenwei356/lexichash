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
	"runtime"
	"sync"

	"github.com/shenwei356/lexichash"
	"github.com/shenwei356/lexichash/tree"
)

// Strands could be used to output strand for a reverse complement flag
var Strands = [2]byte{'+', '-'}

// ErrKConcurrentInsert occurs when calling Insert during calling BatchInsert.
var ErrKConcurrentInsert = errors.New("index: concurrent insertion")

// Index creates LexicHash index for mutitple reference sequences
// and supports searching with a query sequence.
type Index struct {
	lh          *lexichash.LexicHash
	k           uint8
	batchInsert bool

	// each record of the k-mer value is an uint64
	//  ref idx: 26 bits
	//  pos:     36 bits (0-based position)
	//  strand:   2 bits
	Trees []*tree.Tree

	IDs         [][]byte     // IDs of the reference genomes
	RefSeqInfos []RefSeqInfo // Reference sequence basic information
	i           uint32       // curent index, for inserting a new ref seq

	path string // path of the index directory

	KmerLocations [][]uint64 // mask locations of each reference genomes
}

// GenomeInfo is a struct to store some basic information of a ref seq
type RefSeqInfo struct {
	GenomeSize int   // bases of all sequences
	Len        int   // length of contatenated sequences
	NumSeqs    int   // number of sequences
	SeqSizes   []int // sizes of sequences
}

// NewIndex ceates a new Index.
// nMasks better be >= 1024 and better be power of 4,
// i.e., 4, 16, 64, 256, 1024, 4096 ...
func NewIndex(k int, nMasks int) (*Index, error) {
	return NewIndexWithSeed(k, nMasks, 1)
}

// NewIndexWithSeed ceates a new Index with given seed.
// nMasks better be >= 1024 and better be power of 4,
// i.e., 4, 16, 64, 256, 1024, 4096 ...
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
// BatchInsert is recommended, which is faster and safer.
func (idx *Index) Insert(id []byte, s []byte, seqSize int, seqSizes []int) error {
	if idx.batchInsert {
		return ErrKConcurrentInsert
	}

	// compute regions to skip
	var skipRegions [][2]int
	if len(seqSizes) > 1 {
		k := idx.K()
		skipRegions = make([][2]int, len(seqSizes)-1)
		var n int // len of concatenated seqs
		for i, s := range seqSizes {
			if i > 0 {
				skipRegions[i-1] = [2]int{n, n + k - 1}
				n += k - 1
			}
			n += s
		}
	}

	_kmers, locses, err := idx.lh.Mask(s, skipRegions)
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

// PoolRefSeq is the object pool of RefSeq.
var PoolRefSeq = &sync.Pool{New: func() interface{} {
	return &RefSeq{
		ID:         make([]byte, 0, 128),
		Seq:        make([]byte, 0, 10<<20),
		RefSeqSize: 0,
		SeqSizes:   make([]int, 0, 128),
	}
}}

// Reset resets RefSeq.
func (r *RefSeq) Reset() {
	r.ID = r.ID[:0]
	r.Seq = r.Seq[:0]
	r.RefSeqSize = 0
	r.SeqSizes = r.SeqSizes[:0]
}

// MaskResult represents a mask result, it's only used in BatchInsert.
type MaskResult struct {
	ID     []byte
	Kmers  *[]uint64
	Locses *[][]int

	RefSeqSize int   // genome size
	SeqSizes   []int // lengths of each sequences
}

// BatchInsert inserts reference sequences in parallel.
// It returns:
//
//	chan RefSeq, for sending sequence.
//	sync.WaitGroup, for wait all masks being computed.
//	chan int, for waiting all the insertions to be done.
//
// Example:
//
//	input, done := BatchInsert()
//
//	refseq := index.PoolRefSeq.Get().(*index.RefSeq)
//	refseq.Reset()
//
//	// record is a fastx.Record//
//	refseq.ID = append(refseq.ID, record.ID...)
//	refseq.Seq = append(refseq.Seq, record.Seq.Seq...)
//	refseq.SeqSizes = append(refseq.SeqSizes, len(record.Seq.Seq))
//	refseq.RefSeqSize = len(record.Seq.Seq)
//
//	input <- refseq
//
//	close(input)
//	<- done
//
// Multiple sequences can also be concatenated with (K-1) N's for being a single sequence.
// In this case, k-mers around the (K-1) N's regions will be ignored.
func (idx *Index) BatchInsert() (chan *RefSeq, chan int) {
	if idx.batchInsert {
		panic(ErrKConcurrentInsert)
	}
	idx.batchInsert = true

	input := make(chan *RefSeq, Threads)
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
			go func(ref *RefSeq) {
				// compute regions to skip
				var skipRegions [][2]int
				if len(ref.SeqSizes) > 1 {
					k := idx.K()
					skipRegions = make([][2]int, len(ref.SeqSizes)-1)
					var n int // len of concatenated seqs
					for i, s := range ref.SeqSizes {
						if i > 0 {
							skipRegions[i-1] = [2]int{n, n + k - 1}
							n += k - 1
						}
						n += s
					}
				}

				// capture k-mers
				_kmers, locses, err := idx.lh.Mask(ref.Seq, skipRegions)
				if err != nil {
					panic(err)
				}

				m := poolMaskResult.Get().(*MaskResult)
				m.Kmers = _kmers
				m.Locses = locses
				m.ID = ref.ID
				m.RefSeqSize = ref.RefSeqSize
				m.SeqSizes = m.SeqSizes[:0]
				m.SeqSizes = append(m.SeqSizes, ref.SeqSizes...)
				ch <- m

				PoolRefSeq.Put(ref)

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
