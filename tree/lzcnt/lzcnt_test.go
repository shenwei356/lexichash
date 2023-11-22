package lzcnt

import (
	"math/bits"
	"math/rand"
	"testing"
)

func TestLzcnt64(t *testing.T) {
	g := rand.New(rand.NewSource(1))
	n := 10000

	var v uint64
	var e, r int
	for i := 0; i < n; i++ {
		v = g.Uint64()
		e = bits.LeadingZeros64(v)
		r = Lzcnt64(v)
		if i < 10 {
			t.Logf("v: %064b, e: %d, r: %d\n", v, e, r)
		}
		if e != r {
			t.Errorf("Lzcnt64 error: expected=%d, result=%d", e, r)
			return
		}
	}
}
