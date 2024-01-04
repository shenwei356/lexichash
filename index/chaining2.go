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
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.

package index

import (
	"sync"
)

// Chaining2Options contains all options in chaining.
type Chaining2Options struct {
	MaxGap   int
	MinScore int

	// only used in Chain2
	MaxDistance int
	Band        int // only check i in range of  i − A < j < i
}

// DefaultChaining2Options is the defalt vaule of Chaining2Option.
var DefaultChaining2Options = Chaining2Options{
	MaxGap:   32,
	MinScore: 7,

	MaxDistance: 50,
	Band:        20,
}

// Chainer2 is an object for chaining the anchors in two similar sequences.
// Different from Chainer, Chainer2 find chains with no overlaps.
// Anchors/seeds/substrings in Chainer2 is denser than those in Chainer,
// and the chaining score function is also much simpler, only considering
// the lengths of anchors and gaps between them.
type Chainer2 struct {
	options *Chaining2Options

	// scores        []int
	maxscores     []int
	maxscoresIdxs []int
	visited       []bool
	bounds        []int // 4 * chains
}

// NewChainer creates a new chainer.
func NewChainer2(options *Chaining2Options) *Chainer2 {
	c := &Chainer2{
		options: options,

		// scores:        make([]int, 0, 10240),
		maxscores:     make([]int, 0, 10240),
		maxscoresIdxs: make([]int, 0, 10240),
		visited:       make([]bool, 0, 10240),
		bounds:        make([]int, 32),
	}
	return c
}

// RecycleChainingResult reycles the chaining paths.
// Please remember to call this after using the results.
func RecycleChaining2Result(chains *[]*[]int) {
	for _, chain := range *chains {
		poolChain.Put(chain)
	}
	poolChains.Put(chains)
}

var poolChains2 = &sync.Pool{New: func() interface{} {
	tmp := make([]*[]int, 0, 8)
	return &tmp
}}

var poolChain2 = &sync.Pool{New: func() interface{} {
	tmp := make([]int, 0, 32)
	return &tmp
}}

// Chain finds the possible chain paths.
// Please remember to call RecycleChainingResult after using the results.
// Returned results:
//  1. Paths,
//  2. The number of matched bases.
//  3. The number of aligned bases.
func (ce *Chainer2) Chain(subs *[]*SubstrPair) (*[]*[]int, int, int) {
	n := len(*subs)

	if n == 1 { // for one seed, just check the seed weight
		paths := poolChains2.Get().(*[]*[]int)
		*paths = (*paths)[:0]

		sub := (*subs)[0]
		if sub.Len >= ce.options.MinScore { // the length of anchor
			path := poolChain2.Get().(*[]int)
			*path = (*path)[:0]

			*path = append(*path, 0)

			*paths = append(*paths, path)

			return paths, sub.Len, sub.Len
		}

		return paths, 0, 0
	}

	var i, _b, j, k int
	band := ce.options.Band // band size of banded-DP

	// a list for storing score matrix, the size is band * len(seeds pair)
	// scores := ce.scores[:0]
	// size := n * (band + 1)
	// for k = 0; k < size; k++ {
	// 	scores = append(scores, 0)
	// }

	// reused objects

	// the maximum score for each seed, the size is n
	maxscores := ce.maxscores[:0]
	// index of previous seed, the size is n. pointers for backtracking.
	maxscoresIdxs := ce.maxscoresIdxs[:0]
	// for markinng check anchors
	visited := ce.visited[:0]

	// initialize
	maxscores = append(maxscores, (*subs)[0].Len)
	maxscoresIdxs = append(maxscoresIdxs, 0)
	visited = append(visited, false)

	// compute scores
	var s, m, M, d, g int
	var mj, Mi int
	var a, b *SubstrPair
	maxGap := ce.options.MaxGap
	maxDistance := ce.options.MaxDistance
	// scores[0] = (*subs)[0].Len
	for i = 1; i < n; i++ {
		a = (*subs)[i] // current seed/anchor
		k = band * i   // index of current seed in the score matrix

		// just initialize the max score, which comes from the current seed
		m, mj = a.Len, i
		// scores[k] = m

		for _b = 1; _b <= band; _b++ { // check previous $band seeds
			j = i - _b // index of the previous seed
			if j < 0 {
				break
			}

			b = (*subs)[j] // previous seed/anchor
			k++            // index of previous seed in the score matrix

			if b.TBegin > a.TBegin { // filter out messed/crossed anchors
				continue
			}

			d = distance2(a, b)
			if d > maxDistance { // limit the distance. necessary?
				continue
			}

			g = gap2(a, b)
			if g > maxGap { // limit the gap. necessary?
				continue
			}

			s = maxscores[j] + b.Len - g // compute the socre
			// scores[k] = s                // necessary?

			if s >= m { // update the max score of current seed/anchor
				m = s
				mj = j
			}
		}

		maxscores = append(maxscores, m)          // save the max score of the whole
		maxscoresIdxs = append(maxscoresIdxs, mj) // save where the max score comes from
		visited = append(visited, false)          // append marks

		if m > M { // the biggest score in the whole score matrix
			M, Mi = m, i
		}
	}

	// print the score matrix
	// fmt.Printf("i\tpair-i\tiMax\tj:scores\n")
	// for i = 0; i < n; i++ {
	// 	fmt.Printf("%d\t%s\t%d", i, (*subs)[i], maxscoresIdxs[i])
	// 	// k = i * band
	// 	// for _b = 0; _b <= band; _b++ {
	// 	// 	if i-_b >= 0 {
	// 	// 		fmt.Printf("\t%3d:%-4d", i-_b, scores[k])
	// 	// 	}

	// 	// 	k++
	// 	// }
	// 	fmt.Printf("\n")
	// }

	// backtrack

	var firstAnchorOfAChain bool
	minScore := ce.options.MinScore

	paths := poolChains.Get().(*[]*[]int)
	*paths = (*paths)[:0]
	var nMatchedBases, nAlignedBases int

	// ---------------- the chain with the highest score  -----------------

	path := poolChain.Get().(*[]int)
	*path = (*path)[:0]

	i = Mi
	firstAnchorOfAChain = true
	var qb, qe, tb, te int // the bound
	var sub *SubstrPair = (*subs)[i]
	qe, te = sub.QBegin+sub.Len, sub.TBegin+sub.Len // end
	var beginOfNextAnchor int                       // for counting matched bases
	for {
		if firstAnchorOfAChain && maxscores[i] < minScore { // the highest score is not good enough
			return paths, 0, 0 // an empty path
		}

		j = maxscoresIdxs[i] // previous seed

		*path = append(*path, i) // record the seed
		visited[i] = true        // mark as visited
		sub = (*subs)[i]
		if firstAnchorOfAChain {
			firstAnchorOfAChain = false
			nMatchedBases += sub.Len
		} else {
			if sub.QBegin+sub.Len >= beginOfNextAnchor {
				nMatchedBases += beginOfNextAnchor - sub.QBegin
			} else {
				nMatchedBases += sub.Len
			}
		}
		beginOfNextAnchor = sub.QBegin

		if i == j { // the path starts here
			break
		}
		i = j
	}
	reverseInts(*path)
	*paths = append(*paths, path)

	// fmt.Printf("first max: i=%d, score=%d, \n", Mi, M)

	sub = (*subs)[i]
	qb, tb = sub.QBegin, sub.TBegin // begin

	// intervals for checking overlap
	bounds := ce.bounds[:0]
	bounds = append(bounds, qb)
	bounds = append(bounds, qe)
	bounds = append(bounds, tb)
	bounds = append(bounds, te)
	nAlignedBases += qe - qb + 1
	// fmt.Printf("new bound: (%d, %d) vs (%d, %d)\n", qb, qe, tb, te)

	// ------------------------- other chains ----------------------------

	path = poolChain.Get().(*[]int)
	*path = (*path)[:0]

	i = n - 1 // the biggest unvisited i
	firstAnchorOfAChain = true
	var overlapped bool
	var biggestUnvisitedI int
	var computeBiggestUnvisitedI bool
	var nb, bi, bj int // index of bounds
	for {
		// find the next highest score
		M, Mi = 0, i
		computeBiggestUnvisitedI = true
		for ; i >= 0; i-- {
			if visited[i] {
				continue
			}
			if computeBiggestUnvisitedI {
				biggestUnvisitedI = i
				computeBiggestUnvisitedI = false
			}
			m = maxscores[i]
			if m > M {
				M, Mi = m, i
			}
		}
		if M < minScore { // no valid anchors
			break
		}
		i = Mi

		// fmt.Printf("next max: i=%d, score=%d, j=%d, biggestUnvisitedI: %d \n",
		// 	Mi, M, maxscoresIdxs[i], biggestUnvisitedI)

		// find valid anchor for this chain
		for {
			j = maxscoresIdxs[i] // previous seed
			// fmt.Printf("  i:%d, j:%d\n", i, j)
			if visited[j] { // curent seed is abandoned
				visited[i] = true // mark as visited
				// fmt.Printf("  %d is abandoned, j=%d\n", i, j)
				break
			}

			// check if an anchor overlaps with previous chains
			//
			// Query
			// |        te  / (OK)
			// |        |  /
			// |(NO)/   |____qe
			// |   /   /
			// |qb____/    / (NO)
			// |   /  |   /
			// |OK/   |tb
			// o-------------------- Ref
			//
			sub = (*subs)[i]
			overlapped = false
			nb = len(bounds) / 4
			for bi = 0; bi < nb; bi++ {
				bj = bi << 2
				if !((sub.QBegin > bounds[bj+1] && sub.TBegin > bounds[bj+3]) || // top right
					(sub.QBegin+sub.Len < bounds[bj] && sub.TBegin+sub.Len < bounds[bj+2])) { // bottom left
					overlapped = true
					break
				}
			}
			if overlapped {
				// fmt.Printf("  %d (%s) is overlapped previous chain, j=%d\n", i, *sub, j)
				visited[i] = true // mark as visited
				i = j             // check the previous anchor
				continue
			}

			// fmt.Printf("  add %d (%s)\n", i, *sub)
			*path = append(*path, i) // record the seed
			sub = (*subs)[i]
			if firstAnchorOfAChain {
				firstAnchorOfAChain = false

				qe, te = sub.QBegin+sub.Len, sub.TBegin+sub.Len // end
				qb, tb = sub.QBegin, sub.TBegin                 // begin

				nMatchedBases += sub.Len
			} else {
				qb, tb = sub.QBegin, sub.TBegin // begin

				if sub.QBegin+sub.Len >= beginOfNextAnchor {
					nMatchedBases += beginOfNextAnchor - sub.QBegin
				} else {
					nMatchedBases += sub.Len
				}
			}
			beginOfNextAnchor = sub.QBegin

			if i == j { // the path starts here
				visited[i] = true // mark as visited
				reverseInts(*path)
				*paths = append(*paths, path)

				// add the interval of this chain
				bounds = append(bounds, qb)
				bounds = append(bounds, qe)
				bounds = append(bounds, tb)
				bounds = append(bounds, te)
				nAlignedBases += qe - qb + 1
				// fmt.Printf("  new bound: (%d, %d) vs (%d, %d)\n", qb, qe, tb, te)

				path = poolChain.Get().(*[]int)
				*path = (*path)[:0]

				firstAnchorOfAChain = true
				break
			}

			visited[i] = true // mark as visited
			i = j
		}

		i = biggestUnvisitedI
	}

	return paths, nMatchedBases, nAlignedBases
}

func distance2(a, b *SubstrPair) int {
	q := a.QBegin - b.QBegin
	t := a.TBegin - b.TBegin
	if q > t {
		return q
	}
	return t
}

func gap2(a, b *SubstrPair) int {
	g := a.QBegin - b.QBegin - (a.TBegin - b.TBegin)
	if g < 0 {
		return -g
	}
	return g
}
