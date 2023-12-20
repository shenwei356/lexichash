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
)

type ChainingOption struct {
	MaxGap float64
}

func chaining(r *SearchResult, cf *ChainingOption) ([][]int, float64) {
	subs := r.Subs
	n := len(*subs)

	scores := make([]float64, n*(n+1)>>1) // a list for storing triangular score matrix
	maxscores := make([]float64, n)       // the maximum score for each seed
	maxscoresIdxs := make([]int, n)       // index of previous seed

	var i, j, j0, k, mj int
	// initialize
	for i, b := range *subs { // j == i, means a path starting from this seed
		j0 = i * (i + 1) / 2
		k = j0 + i
		scores[k] = seedWeight(b)
	}
	maxscores[0] = scores[0]
	maxscoresIdxs[0] = 0

	// compute scores
	var s, m, d float64
	var a, b *SubstrPair
	for i = 1; i < n; i++ {
		j0 = i * (i + 1) / 2

		// just initialize the max score, which comes from the curent seed
		k = j0 + i // starting with seed i
		m = scores[k]
		mj = i

		for j = 0; j < i; j++ { // try all previous seeds
			k = j0 + j
			a, b = (*subs)[i], (*subs)[j]

			d = distance(a, b)
			if d > cf.MaxGap {
				continue
			}

			s = maxscores[i-1] + seedWeight(b) - distanceScore(d) - gapScore(a, b)
			scores[k] = s

			if s > m { // update the max score
				m = s
				mj = j
			}
		}
		maxscores[i] = m // save the max score
		maxscoresIdxs[i] = mj
	}
	// print the score matrix
	// for i = 0; i < n; i++ {
	// 	fmt.Printf("%d", i+1)
	// 	for j = 0; j <= i; j++ {
	// 		k = i*(i+1)/2 + j
	// 		fmt.Printf("\t%5.1f", scores[k])
	// 	}
	// 	fmt.Printf("\n")
	// }

	// backtrack
	visited := make([]bool, n)
	paths := make([][]int, 0, n)
	var first bool
	var sumMaxScore float64

	path := make([]int, 0, n)
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

		path = append(path, i) // record the seed
		visited[i] = true      // mark as visited
		j = maxscoresIdxs[i]   // previous seed
		if first {
			sumMaxScore += maxscores[i]
			first = false
		}
		if i != j {
			i = j
		} else { // the path starts here
			reverseInts(path)
			paths = append(paths, path)

			path = make([]int, 0, n)
			i = n - 1 // re-track from the end
			first = true
		}
	}
	// fmt.Println("paths:", paths, sumMaxScore)
	return paths, sumMaxScore
}

func seedWeight(b *SubstrPair) float64 {
	return 0.1 * float64(b.Len) * float64(b.Len)
}

func distance(a, b *SubstrPair) float64 {
	return math.Max(math.Abs(float64(a.QBegin-b.QBegin)), math.Abs(float64(a.TBegin-b.TBegin)))
}

func distanceScore(d float64) float64 {
	d = 0.01 * d
	if d > 10 {
		return 10
	}
	return d
}

func gapScore(a, b *SubstrPair) float64 {
	gap := math.Abs(float64(a.QBegin - b.QBegin - a.TBegin + b.TBegin))
	if gap == 0 {
		return 0
	}
	return 0.01*gap + 0.5*math.Log2(gap)
}

func reverseInts(s []int) {
	for i, j := 0, len(s)-1; i < j; i, j = i+1, j-1 {
		s[i], s[j] = s[j], s[i]
	}
}
