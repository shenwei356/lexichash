//go:build !amd64 && !386

package lzcnt

import "math/bits"

var lzcntFuncs = []lzcntImpl{
	{bits.LeadingZeros64, "generic", true},
}
