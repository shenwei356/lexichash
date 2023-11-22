//go:build ignore
// +build ignore

package main

import (
	. "github.com/mmcloughlin/avo/build"
)

func main() {
	TEXT("lzcnt64", NOSPLIT, "func(v uint64) int")
	Doc("Using lzcnt to count Leading zeros of an uint64.")
	v := Load(Param("v"), GP64())
	n := GP64()
	LZCNTQ(v, n)
	Store(n, ReturnIndex(0))
	RET()
	Generate()
}
