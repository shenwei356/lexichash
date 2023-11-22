package lzcnt

type lzcnt64Impl struct {
	function  func(v uint64) int
	name      string
	available bool
}

var Lzcnt64 = func() func(v uint64) int {
	for _, f := range lzcntFuncs {
		if f.available {
			return f.function
		}
	}

	panic("no implementation available")
}()
