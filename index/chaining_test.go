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
	"testing"
)

func TestChaining(t *testing.T) {
	subs := []*SubstrPair{
		// {QBegin: 419, TBegin: 1419, Len: 31},
		// {QBegin: 450, TBegin: 3638447, Len: 31},
		// {QBegin: 547, TBegin: 3638544, Len: 31},
		{QBegin: 10, TBegin: 1314, Len: 20},

		{QBegin: 60, TBegin: 3395374, Len: 15},
		{QBegin: 70, TBegin: 3395384, Len: 15},

		{QBegin: 50, TBegin: 950, Len: 31},
		{QBegin: 79, TBegin: 3637976, Len: 31},
		{QBegin: 100, TBegin: 3637997, Len: 31},
		{QBegin: 519, TBegin: 1419, Len: 31},
		{QBegin: 550, TBegin: 3638447, Len: 31},
		{QBegin: 647, TBegin: 3638544, Len: 31},

		{QBegin: 111, TBegin: 1146311, Len: 31},
		{QBegin: 136, TBegin: 1146336, Len: 31},
		{QBegin: 138, TBegin: 1146338, Len: 31},
		{QBegin: 139, TBegin: 1146339, Len: 31},
		{QBegin: 264, TBegin: 1146464, Len: 31},
		{QBegin: 1479, TBegin: 1147679, Len: 31},
		{QBegin: 1484, TBegin: 1147684, Len: 31},
		{QBegin: 1543, TBegin: 1147743, Len: 31},
		{QBegin: 1566, TBegin: 1147766, Len: 31},
		{QBegin: 1919, TBegin: 1148119, Len: 31},
	}
	tmp := []*SearchResult{
		{
			IdIdx: 0,
			Subs:  &subs,
		},
	}
	rs := &tmp

	cf := &DefaultChainingOptions

	chainer := NewChainer(cf)
	for _, r := range *rs {
		paths, sumMaxScore := chainer.Chain(r)

		t.Logf("sum score: %f, paths:\n", sumMaxScore)
		for _, p := range *paths {
			t.Logf("  %d\n", *p)
		}

		RecycleChainingResult(paths)
	}
}
