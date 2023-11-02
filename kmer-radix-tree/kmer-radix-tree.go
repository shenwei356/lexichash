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

package tree

import (
	"sort"
)

// leafNode is used to represent a value
type leafNode struct {
	key uint64 // bases in the node
	k   uint8  // lenght of bases of the key
	val []uint64
}

// edge is used to represent an edge node
type edge struct {
	label uint8 // a base in 2-bit code
	node  *node
}

type node struct {
	leaf *leafNode // optional

	prefix uint64 // prefix of current node
	k      uint8  // lenght of bases of the prefix

	edges edges // need to be sorted for fast query
}

func (n *node) isLeaf() bool {
	return n.leaf != nil
}

func (n *node) addEdge(e edge) {
	num := len(n.edges)
	idx := sort.Search(num, func(i int) bool { // find the position to insert
		return n.edges[i].label >= e.label
	})
	n.edges = append(n.edges, edge{})
	copy(n.edges[idx+1:], n.edges[idx:])
	n.edges[idx] = e
}

func (n *node) updateEdge(label uint8, node *node) {
	num := len(n.edges)
	idx := sort.Search(num, func(i int) bool {
		return n.edges[i].label >= label
	})
	if idx < num && n.edges[idx].label == label {
		n.edges[idx].node = node
		return
	}
	panic("replacing missing edge")
}

func (n *node) getEdge(label uint8) *node {
	num := len(n.edges)
	idx := sort.Search(num, func(i int) bool {
		return n.edges[i].label >= label
	})
	if idx < num && n.edges[idx].label == label {
		return n.edges[idx].node
	}
	return nil
}

func (n *node) delEdge(label uint8) {
	num := len(n.edges)
	idx := sort.Search(num, func(i int) bool {
		return n.edges[i].label >= label
	})
	if idx < num && n.edges[idx].label == label {
		copy(n.edges[idx:], n.edges[idx+1:])
		n.edges[len(n.edges)-1] = edge{}
		n.edges = n.edges[:len(n.edges)-1]
	}
}

type edges []edge

func (e edges) Len() int {
	return len(e)
}

func (e edges) Less(i, j int) bool {
	return e[i].label < e[j].label
}

func (e edges) Swap(i, j int) {
	e[i], e[j] = e[j], e[i]
}

func (e edges) Sort() {
	sort.Sort(e)
}

type Tree struct {
	root *node
	size int
}

// Tree implements a radix tree for k-mer query
func New() *Tree {
	t := &Tree{root: &node{}}
	return t
}

// Len returns the number of elements
func (t *Tree) Len() int {
	return t.size
}

// Insert is used to add a newentry or update
// an existing entry. Returns true if an existing record is updated.
func (t *Tree) Insert(key uint64, k uint8, v uint64) ([]uint64, bool) {
	var parent *node
	n := t.root
	key0 := key
	k0 := k
	search := key
	for {
		// Handle key exhaustion
		if k == 0 {
			if n.isLeaf() {
				old := n.leaf.val
				if n.leaf.val == nil {
					n.leaf.val = make([]uint64, 0, 1)
				}
				n.leaf.val = append(n.leaf.val, v)
				return old, true
			}

			n.leaf = &leafNode{
				key: key0,
				k:   k0,
				val: []uint64{v},
			}
			t.size++
			return nil, false
		}

		// Look for the edge
		parent = n
		n = n.getEdge(KmerBaseAt(search, k, 0))

		// No edge, create one
		if n == nil {
			e := edge{
				label: KmerBaseAt(search, k, 0),
				node: &node{
					leaf: &leafNode{
						key: key0,
						k:   k0,
						val: []uint64{v},
					},
					prefix: search,
					k:      k,
				},
			}
			parent.addEdge(e)
			t.size++
			return nil, false
		}

		// Determine longest prefix of the search key on match
		commonPrefix := KmerLongestPrefix(search, n.prefix, k, n.k)
		if commonPrefix == n.k {
			search = KmerSuffix(search, k, commonPrefix) // left bases
			k = k - commonPrefix                         // need to update it
			continue
		}

		// Split the node
		t.size++
		child := &node{
			prefix: KmerPrefix(search, k, commonPrefix), // prefix
			k:      commonPrefix,
		}
		parent.updateEdge(KmerBaseAt(search, k, 0), child)

		// Restore the existing node
		child.addEdge(edge{
			label: KmerBaseAt(n.prefix, n.k, commonPrefix),
			node:  n,
		})
		n.prefix = KmerSuffix(n.prefix, n.k, commonPrefix)
		n.k = n.k - commonPrefix

		// Create a new leaf node
		leaf := &leafNode{
			key: key0,
			k:   k0,
			val: []uint64{v},
		}

		// If the new key is a subset, add to this node
		search = KmerSuffix(search, k, commonPrefix)
		k = k - commonPrefix
		if k == 0 {
			child.leaf = leaf
			return nil, false
		}

		// Create a new edge for the node
		child.addEdge(edge{
			label: KmerBaseAt(search, k, 0),
			node: &node{
				leaf:   leaf,
				prefix: search,
				k:      k,
			},
		})
		return nil, false
	}
}

// Get is used to lookup a specific key, returning
// the value and if it was found
func (t *Tree) Get(key uint64, k uint8) ([]uint64, bool) {
	n := t.root
	search := key
	for {
		// Check for key exhaution
		if k == 0 {
			if n.isLeaf() {
				return n.leaf.val, true
			}
			break
		}

		// Look for an edge
		n = n.getEdge(KmerBaseAt(search, k, 0))
		if n == nil {
			break
		}

		// Consume the search prefix
		if KmerHasPrefix(search, n.prefix, k, n.k) {
			search = KmerSuffix(search, k, n.k)
			k = k - n.k
		} else {
			break
		}
	}
	return nil, false
}

// LongestPrefix is like Get, but instead of an
// exact match, it will return the longest prefix match.
func (t *Tree) LongestPrefix(key uint64, k uint8) (uint64, uint8, []uint64, bool) {
	var last *leafNode
	n := t.root
	search := key
	for {
		// Look for a leaf node
		if n.isLeaf() {
			last = n.leaf
		}

		// Check for key exhaution
		if k == 0 {
			break
		}

		// Look for an edge
		n = n.getEdge(KmerBaseAt(search, k, 0))
		if n == nil {
			break
		}

		// Consume the search prefix
		if KmerHasPrefix(search, n.prefix, k, n.k) {
			search = KmerSuffix(search, k, n.k)
			k = k - n.k
		} else {
			break
		}
	}
	if last != nil {
		return last.key, last.k, last.val, true
	}
	return 0, 0, nil, false
}

// WalkFn is used when walking the tree. Takes a
// key and value, returning if iteration should
// be terminated.
type WalkFn func(key uint64, k uint8, v []uint64) bool

// Walk is used to walk the tree
func (t *Tree) Walk(fn WalkFn) {
	recursiveWalk(t.root, fn)
}

// recursiveWalk is used to do a pre-order walk of a node
// recursively. Returns true if the walk should be aborted
func recursiveWalk(n *node, fn WalkFn) bool {
	// Visit the leaf values if any
	if n.leaf != nil && fn(n.leaf.key, n.leaf.k, n.leaf.val) {
		return true
	}

	// Recurse on the children
	i := 0
	k := len(n.edges) // keeps track of number of edges in previous iteration
	for i < k {
		e := n.edges[i]
		if recursiveWalk(e.node, fn) {
			return true
		}
		// It is a possibility that the WalkFn modified the node we are
		// iterating on. If there are no more edges, mergeChild happened,
		// so the last edge became the current node n, on which we'll
		// iterate one last time.
		if len(n.edges) == 0 {
			return recursiveWalk(n, fn)
		}
		// If there are now less edges than in the previous iteration,
		// then do not increment the loop index, since the current index
		// points to a new edge. Otherwise, get to the next index.
		if len(n.edges) >= k {
			i++
		}
		k = len(n.edges)
	}
	return false
}
