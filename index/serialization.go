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

package index

import (
	"encoding/binary"
	"errors"
)

var be = binary.BigEndian

var Magic = [8]byte{'l', 'e', 'x', 'i', 'c', 'i', 'd', 'x'}

var MainVersion uint8 = 0
var MinorVersion uint8 = 1

// ErrInvalidFileFormat means invalid file format.
var ErrInvalidFileFormat = errors.New("lexichash index: invalid binary format")

// ErrBrokenFile means the file is not complete.
var ErrBrokenFile = errors.New("lexichash index: broken file")

// ErrVersionMismatch means version mismatch between files and program
var ErrVersionMismatch = errors.New("lexichash: version mismatch")

// NewFromFile creates a LexicHash from a file
func NewFromPath(file string) (*Index, error) {

	return nil, nil
}

// WriteToFile writes a LexicHash to a file, optional with file extension of .gz, .xz, .zst, .bz2.
func (idx *Index) WriteToPath(file string) (int, error) {

	return 0, nil
}
