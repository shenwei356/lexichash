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

package lexichash

import (
	"encoding/binary"
	"errors"
	"io"
	"sync"

	"github.com/shenwei356/xopen"
)

var be = binary.BigEndian

var Magic = [8]byte{'k', 'm', 'e', 'r', 'm', 'a', 's', 'k'}

var MainVersion uint8 = 0
var MinorVersion uint8 = 1

// ErrInvalidFileFormat means invalid file format.
var ErrInvalidFileFormat = errors.New("lexichash: invalid binary format")

// ErrBrokenFile means the file is not complete.
var ErrBrokenFile = errors.New("lexichash: broken file")

// ErrVersionMismatch means version mismatch between files and program.
var ErrVersionMismatch = errors.New("lexichash: version mismatch")

// NewFromFile creates a LexicHash from a file.
func NewFromFile(file string) (*LexicHash, error) {
	fh, err := xopen.Ropen(file)
	if err != nil {
		return nil, err
	}
	defer fh.Close()

	return Read(fh)
}

// WriteToFile writes a LexicHash to a file,
// optional with file extensions of .gz, .xz, .zst, .bz2.
func (lh *LexicHash) WriteToFile(file string) (int, error) {
	outfh, err := xopen.Wopen(file)
	if err != nil {
		return 0, err
	}
	defer outfh.Close()

	return lh.Write(outfh)
}

// Write writes a LexicHash.
//
// Header (32 bytes):
//
//	Magic number, 8 bytes, kmermask
//	Main and minor versions, 2 bytes
//	K, 1 byte
//	Blank, 5 bytes
//	Seed: 8 bytes
//	Number of masks: 8 bytes
//
// Data: k-mers.
//
//	K-mers in uint64, 8*$(the number of maskes)
func (lh *LexicHash) Write(w io.Writer) (int, error) {
	var N int // the number of bytes.
	var err error

	// 8-byte magic number
	err = binary.Write(w, be, Magic)
	if err != nil {
		return N, err
	}
	N += 8

	// 8-byte meta info
	err = binary.Write(w, be, [8]uint8{MainVersion, MinorVersion, uint8(lh.K)})
	if err != nil {
		return N, err
	}
	N += 8

	// 8-byte seed
	err = binary.Write(w, be, uint64(lh.Seed))
	if err != nil {
		return N, err
	}
	N += 8

	// 8-byte the number of masks
	err = binary.Write(w, be, uint64(len(lh.Masks)))
	if err != nil {
		return N, err
	}
	N += 8

	data := make([]byte, 8*len(lh.Masks))
	var i int
	for _, mask := range lh.Masks {
		be.PutUint64(data[i:i+8], mask)
		i += 8
	}
	_, err = w.Write(data)
	if err != nil {
		return N, err
	}
	N += len(data)

	return N, nil
}

// Read reads a LexiHash from an io.Reader.
func Read(r io.Reader) (*LexicHash, error) {
	buf := make([]byte, 64)

	var err error
	var n int

	// check the magic number
	n, err = io.ReadFull(r, buf[:8])
	if err != nil {
		return nil, err
	}
	if n < 8 {
		return nil, ErrBrokenFile
	}
	same := true
	for i := 0; i < 8; i++ {
		if Magic[i] != buf[i] {
			same = false
			break
		}
	}
	if !same {
		return nil, ErrInvalidFileFormat
	}

	// read metadata
	n, err = io.ReadFull(r, buf[:8])
	if err != nil {
		return nil, err
	}
	if n < 8 {
		return nil, ErrBrokenFile
	}
	// check compatibility
	if MainVersion != buf[0] {
		return nil, ErrVersionMismatch
	}
	// check k-mer size
	if buf[2] > 32 {
		return nil, ErrKOverflow
	}

	lh := &LexicHash{K: int(buf[2])}

	// the seed
	_, err = io.ReadFull(r, buf[:8])
	if err != nil {
		return nil, err
	}
	lh.Seed = int64(be.Uint64(buf[:8]))

	// the number of k-mers
	_, err = io.ReadFull(r, buf[:8])
	if err != nil {
		return nil, err
	}
	nKmers := be.Uint64(buf[:8])
	masks := make([]uint64, nKmers)

	buf8 := make([]byte, 8)
	var nReaded int
	for i := 0; i < int(nKmers); i++ {
		nReaded, err = io.ReadFull(r, buf8)
		if err != nil {
			return nil, err
		}
		if nReaded < 8 {
			return nil, ErrBrokenFile
		}

		masks[i] = be.Uint64(buf8)
	}
	lh.Masks = masks

	// do not forgot the buffer

	lh.poolKmers = &sync.Pool{New: func() interface{} {
		kmers := make([]uint64, len(masks))
		return &kmers
	}}
	lh.poolLocses = &sync.Pool{New: func() interface{} {
		locses := make([][]int, len(masks))
		for i := range locses {
			locses[i] = make([]int, 1)
		}
		return &locses
	}}
	lh.poolHashes = &sync.Pool{New: func() interface{} {
		hashes := make([]uint64, len(masks))
		return &hashes
	}}

	return lh, nil
}
