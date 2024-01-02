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
	"fmt"
	"sync"
)

// ChainingOptions contains all options in chaining.
type Chaining2Options struct {
	MaxGap   int
	MinScore int

	MinDistance int // minimum distance between two seeds

	// only used in Chain2
	MaxDistance int
	Band        int // only check i in range of  i − A < j < i
}

// DefaultChainingOptions is the defalt vaule of ChainingOption.
var DefaultChaining2Options = ChainingOptions{
	MaxGap:   5000,
	MinScore: 40,

	MinDistance: 5,
}

// Chainer is an object for chaining the seeds.
// Some variables like the score tables are re-used,
// so it could help to reduce GC load.
type Chainer2 struct {
	options *Chaining2Options

	scores        []int
	maxscores     []int
	maxscoresIdxs []int
	visited       []bool
}

// NewChainer creates a new chainer.
func NewChainer2(options *Chaining2Options) *Chainer2 {
	c := &Chainer2{
		options: options,

		scores:        make([]int, 0, 1024),
		maxscores:     make([]int, 0, 1024),
		maxscoresIdxs: make([]int, 0, 128),
		visited:       make([]bool, 0, 128),
	}
	return c
}

// RecycleChainingResult reycles the chaining results.
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

// Chain finds the possible seed paths.
// Please remember to call RecycleChainingResult after using the results.
func (ce *Chainer2) Chain(subs *[]*SubstrPair) (*[]*[]int, int) {
	n := len(*subs)

	var sumMaxScore int

	if n == 1 { // for one seed, just check the seed weight
		paths := poolChains2.Get().(*[]*[]int)
		*paths = (*paths)[:0]

		w := (*subs)[0].Len
		if w >= ce.options.MinScore {
			path := poolChain2.Get().(*[]int)
			*path = (*path)[:0]

			*path = append(*path, 0)

			*paths = append(*paths, path)
		}

		return paths, w
	}

	var i, _b, j, k int
	band := ce.options.Band

	// a list for storing score matrix, the size is band * len(seeds pair)
	scores := ce.scores[:0]
	size := n * (band + 1)
	for k = 0; k < size; k++ {
		scores = append(scores, 0)
	}
	// the maximum score for each seed, the size is n
	maxscores := ce.maxscores[:0]
	for i = 0; i < n; i++ {
		maxscores = append(maxscores, 0)
	}
	// index of previous seed, the size is n
	maxscoresIdxs := ce.maxscoresIdxs[:0]
	for i = 0; i < n; i++ {
		maxscoresIdxs = append(maxscoresIdxs, 0)
	}

	// initialize
	maxscores[0] = (*subs)[0].Len
	maxscoresIdxs[0] = 0

	// compute scores
	var s, m, M, d, g int
	var mj, Mi int
	var a, b *SubstrPair
	maxGap := ce.options.MaxGap
	maxDistance := ce.options.MaxDistance
	for i = 1; i < n; i++ {
		a = (*subs)[i]
		k = band * i
		// just initialize the max score, which comes from the current seed
		m = a.Len
		mj = i
		scores[k] = m
		for _b = 1; _b <= band; _b++ { // check previous $band seeds
			j = i - _b
			if j < 0 {
				break
			}
			k++

			b = (*subs)[j]

			if b.TBegin > a.TBegin {
				continue
			}

			d = distance2(a, b)
			if d > maxDistance {
				continue
			}

			g = gap2(a, b)
			if g > maxGap {
				continue
			}

			s = maxscores[j] + b.Len - g
			scores[k] = s

			if s >= m { // update the max score
				m = s
				mj = j
			}
		}
		maxscores[i] = m // save the max score
		maxscoresIdxs[i] = mj

		if m > M { // the biggest score
			M = m
			Mi = i
		}
	}

	// print the score matrix
	fmt.Printf("i\tpair-i\tiMax\tscores\n")
	for i = 0; i < n; i++ {
		fmt.Printf("%d\t%s\t%d", i, (*subs)[i], maxscoresIdxs[i])
		k = i * band
		for _b = 0; _b <= band; _b++ {
			if i-_b >= 0 {
				fmt.Printf("\t%3d:%-4d", i-_b, scores[k])
			}

			k++
		}
		fmt.Printf("\n")
	}

	// backtrack
	visited := ce.visited[:0]
	for i = 0; i < n; i++ {
		visited = append(visited, false)
	}

	var first bool
	minScore := ce.options.MinScore

	paths := poolChains.Get().(*[]*[]int)
	*paths = (*paths)[:0]

	// ---------------- the longest chain -----------------
	path := poolChain.Get().(*[]int)
	*path = (*path)[:0]

	i = Mi
	first = true
	var qb, qe, tb, te int // the bound of the longest chain
	var sub *SubstrPair = (*subs)[i]
	qe, te = sub.QBegin+sub.Len, sub.TBegin+sub.Len
	for {
		if first && maxscores[i] < minScore { // the highest score is not good enough
			return paths, maxscores[i]
		}

		j = maxscoresIdxs[i] // previous seed

		*path = append(*path, i) // record the seed
		visited[i] = true        // mark as visited
		if first {
			sumMaxScore += maxscores[i]
			first = false
		}

		if i == j { // the path starts here
			break
		}
		i = j
	}
	reverseInts(*path)
	*paths = append(*paths, path)

	sub = (*subs)[i]
	qb, tb = sub.QBegin, sub.TBegin

	fmt.Println(qb, qe, tb, te)

	// ---------------- other chains -----------------
	path = poolChain.Get().(*[]int)
	*path = (*path)[:0]

	i = n - 1 // the bigest unvisited i
	var biggestUnvisitedI int
	first = true
	for {
		// find the larget unvisited i
		for ; i >= 0; i-- {
			if !visited[i] {
				biggestUnvisitedI = i
				break
			}
		}
		// all are visited
		if i == -1 {
			break
		}

		if first && maxscores[i] < minScore {
			visited[i] = true
			i--
			continue
		}

		sub = (*subs)[i]
		// Query
		// |        te  / (OK)
		// |        |  /
		// |(NO)/   |____qe
		// |   /   /
		// |qb____/    / (NO)
		// |   /  |   /
		// |OK/   |tb
		// o-------------------- Ref
		if !((sub.QBegin > qe && sub.TBegin > te) || // top right
			(sub.QBegin+sub.Len < qb && sub.TBegin+sub.Len < tb)) { // bottom left
			visited[i] = true
			i--
			continue
		}

		j = maxscoresIdxs[i] // previous seed
		if visited[j] {      // curent seed is abandoned
			i--
			continue
		}

		*path = append(*path, i) // record the seed
		visited[i] = true        // mark as visited
		if first {
			sumMaxScore += maxscores[i]
			first = false
		}
		if i != j {
			i = j
		} else { // the path starts here
			reverseInts(*path)
			*paths = append(*paths, path)

			path = poolChain.Get().(*[]int)
			*path = (*path)[:0]

			i = biggestUnvisitedI - 1
			first = true
		}
	}

	return paths, sumMaxScore
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
