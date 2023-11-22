//go:build amd64 || 386

package lzcnt

import (
	"math/bits"

	"golang.org/x/sys/cpu"
)

var lzcntFuncs = []lzcnt64Impl{
	{lzcnt64, "SSE42", cpu.X86.HasSSE42},
	{bits.LeadingZeros64, "generic", true},
}
