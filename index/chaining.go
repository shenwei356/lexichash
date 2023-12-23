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
	"math"
	"sync"
)

// ChainingOptions contains all options in chaining.
type ChainingOptions struct {
	MaxGap   float64
	MinScore float64
}

// DefaultChainingOptions is the defalt vaule of ChainingOption.
var DefaultChainingOptions = ChainingOptions{
	MaxGap:   5000,
	MinScore: 40,
}

// Chainer is an object for chaining the seeds.
// Some variables like the score tables are re-used,
// so it could help to reduce GC load.
type Chainer struct {
	options *ChainingOptions

	scores        []float64
	maxscores     []float64
	maxscoresIdxs []int
	visited       []bool
}

// NewChainer creates a new chainer.
func NewChainer(options *ChainingOptions) *Chainer {
	c := &Chainer{
		options: options,

		scores:        make([]float64, 0, 128),
		maxscores:     make([]float64, 0, 128),
		maxscoresIdxs: make([]int, 0, 128),
		visited:       make([]bool, 0, 128),
	}
	return c
}

// RecycleChainingResult reycles the chaining results.
// Please remember to call this after using the results.
func RecycleChainingResult(paths *[]*[]int) {
	for _, path := range *paths {
		poolChain.Put(path)
	}
	poolChains.Put(paths)
}

var poolChains = &sync.Pool{New: func() interface{} {
	tmp := make([]*[]int, 0, 8)
	return &tmp
}}

var poolChain = &sync.Pool{New: func() interface{} {
	tmp := make([]int, 0, 32)
	return &tmp
}}

// Chain find the possible seed paths.
// Please remember to call RecycleChainingResult after using the results.
func (ce *Chainer) Chain(r *SearchResult) (*[]*[]int, float64) {
	subs := r.Subs
	n := len(*subs)

	var sumMaxScore float64

	if n == 1 { // for one seed, just check the seed weight
		paths := poolChains.Get().(*[]*[]int)
		*paths = (*paths)[:0]

		w := seedWeight(float64((*subs)[0].Len))
		if w >= ce.options.MinScore {
			path := poolChain.Get().(*[]int)
			*path = (*path)[:0]

			*path = append(*path, 0)

			*paths = append(*paths, path)
		}

		return paths, w
	}

	var i, j, j0, k, mj int

	// a list for storing triangular score matrix, the size is n*(n+1)>>1
	scores := ce.scores[:0]
	for k = 0; k < n*(n+1)>>1; k++ {
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
	for i, b := range *subs { // j == i, means a path starting from this seed
		j0 = i * (i + 1) / 2
		k = j0 + i
		scores[k] = seedWeight(float64(b.Len))
	}
	maxscores[0] = scores[0]
	maxscoresIdxs[0] = 0

	// compute scores
	var s, m, d float64
	var a, b *SubstrPair
	maxGap := ce.options.MaxGap
	for i = 1; i < n; i++ {
		j0 = i * (i + 1) / 2

		// just initialize the max score, which comes from the current seed
		k = j0 + i // starting with seed i
		m = scores[k]
		mj = i

		for j = 0; j < i; j++ { // try all previous seeds
			k = j0 + j
			a, b = (*subs)[i], (*subs)[j]

			d = distance(a, b)
			if d > maxGap {
				continue
			}

			s = maxscores[i-1] + seedWeight(float64(b.Len)) - distanceScore(d) - gapScore(a, b)
			scores[k] = s

			if s >= m { // update the max score
				m = s
				mj = j
			}
		}
		maxscores[i] = m // save the max score
		maxscoresIdxs[i] = mj
	}
	// print the score matrix
	// fmt.Printf("i\tiMax\tscores\n")
	// for i = 0; i < n; i++ {
	// 	fmt.Printf("%d\t%d", i, maxscoresIdxs[i])
	// 	for j = 0; j <= i; j++ {
	// 		k = i*(i+1)/2 + j
	// 		fmt.Printf("\t%6.2f", scores[k])
	// 	}
	// 	fmt.Printf("\n")
	// }

	// backtrack
	visited := ce.visited[:0]
	for i = 0; i < n; i++ {
		visited = append(visited, false)
	}
	paths := poolChains.Get().(*[]*[]int)
	*paths = (*paths)[:0]
	var first bool
	minScore := ce.options.MinScore

	path := poolChain.Get().(*[]int)
	*path = (*path)[:0]
	i = n - 1
	first = true
	for {
		// find the larget unvisited i
		for ; i >= 0; i-- {
			if !visited[i] {
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
			i = n - 1 // re-track from the end
			first = true
		}
	}

	return paths, sumMaxScore
}

func seedWeight(l float64) float64 {
	return 0.1 * l * l
}

func distance(a, b *SubstrPair) float64 {
	return math.Max(math.Abs(float64(a.QBegin-b.QBegin)), math.Abs(float64(a.TBegin-b.TBegin)))
}

func distanceScore(d float64) float64 {
	return 0.01 * d
}

func gapScore(a, b *SubstrPair) float64 {
	gap := math.Abs(float64(a.QBegin - b.QBegin - a.TBegin + b.TBegin))
	if gap == 0 {
		return 0
	}
	return 0.1*gap + 0.5*math.Log2(gap)
}

func reverseInts(s []int) {
	for i, j := 0, len(s)-1; i < j; i, j = i+1, j-1 {
		s[i], s[j] = s[j], s[i]
	}
}
