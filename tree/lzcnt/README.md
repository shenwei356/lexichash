## Generate Go assembly code

Generate Go assembly code with [avo](https://github.com/mmcloughlin/avo),
and perform tests.

```
go run asm-lzcnt64.go -out lzcnt64_amd64.s -stubs lzcnt64.go

go test . -count=1

```
