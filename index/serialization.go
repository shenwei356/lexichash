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
	"bufio"
	"errors"
	"fmt"
	"io/fs"
	"os"
	"path/filepath"
	"runtime"
	"strconv"
	"sync"

	"github.com/shenwei356/lexichash"
	"github.com/shenwei356/lexichash/tree"
	"github.com/shenwei356/util/pathutil"
	"github.com/shenwei356/xopen"
)

var Magic = [8]byte{'l', 'e', 'x', 'i', 'c', 'i', 'd', 'x'}

var MainVersion uint8 = 0
var MinorVersion uint8 = 1

// ErrInvalidFileFormat means invalid file format.
var ErrInvalidFileFormat = errors.New("lexichash index: invalid binary format")

// ErrBrokenFile means the file is not complete.
var ErrBrokenFile = errors.New("lexichash index: broken file")

// ErrVersionMismatch means version mismatch between files and program.
var ErrVersionMismatch = errors.New("lexichash index: version mismatch")

// ErrDirNotEmpty means the output directory is not empty.
var ErrDirNotEmpty = errors.New("lexichash index: output directory not empty")

// ErrDirNotEmpty means the output directory is not empty.
var ErrPWDAsOutDir = errors.New("lexichash index: current directory cant't be the output dir")

// ErrInvalidIndexDir means the path is not a valid index directory.
var ErrInvalidIndexDir = errors.New("lexichash index: invalid index directory")

// ErrTreeFileMissing means some tree files are missing.
var ErrTreeFileMissing = errors.New("lexichash index: some tree files missing")

// IDListFile defines the name of the ID list file
const IDListFile = "IDs.txt"

// MaskFile is the name of the mask file
const MaskFile = "masks.bin"

// TreeDir is the name of director of trees
const TreeDir = "trees"

// TreeFileExt is the file extension of the tree file.
// It would be ".gz", ".xz", ".zst" or ".bz2",
// but they are not recommended when saving a lot of trees,
// as it would assume a lot of RAM.
const TreeFileExt = ".bin"

// WriteToPath writes an index to a directory.
//
// Files:
//
//	Mask file, binary
//	ID list file, plain text
//	Trees directory, binary. Files numbers: 1-5000
func (idx *Index) WriteToPath(outDir string, overwrite bool, threads int) error {
	pwd, _ := os.Getwd()
	if outDir != "./" && outDir != "." && pwd != filepath.Clean(outDir) {
		existed, err := pathutil.DirExists(outDir)
		if err != nil {
			return err
		}
		if existed {
			empty, err := pathutil.IsEmpty(outDir)
			if err != nil {
				return err
			}

			if !empty {
				if overwrite {
					err = os.RemoveAll(outDir)
					if err != nil {
						return err
					}
				} else {
					return ErrDirNotEmpty
				}
			} else {
				err = os.RemoveAll(outDir)
				if err != nil {
					return err
				}
			}
		}
		err = os.MkdirAll(outDir, 0777)
		if err != nil {
			return err
		}
	} else {
		return ErrPWDAsOutDir
	}

	// ID list file
	err := idx.writeIDlist(filepath.Join(outDir, IDListFile))
	if err != nil {
		return err
	}

	// Mask file
	_, err = idx.lh.WriteToFile(filepath.Join(outDir, MaskFile))
	if err != nil {
		return err
	}

	if threads <= 0 {
		threads = runtime.NumCPU()
	}

	// Trees
	var wg sync.WaitGroup
	tokens := make(chan int, threads)

	var idStr, subDir, file string
	for i, t := range idx.Trees {
		idStr = fmt.Sprintf("%04d", i)
		subDir = idStr[len(idStr)-2:]
		file = filepath.Join(outDir, TreeDir, subDir, idStr+TreeFileExt)

		wg.Add(1)
		tokens <- 1

		go func(t *tree.Tree, file string) {
			defer func() {
				wg.Done()
				<-tokens
			}()
			t.WriteToFile(file)
			_, _err := t.WriteToFile(file)
			if _err != nil {
				err = _err
			}
		}(t, file)
	}
	wg.Wait()

	return err
}

// NewFromPath reads an index from a directory.
func NewFromPath(outDir string, threads int) (*Index, error) {
	// ------------- checking directory structure -----------

	ok, err := pathutil.DirExists(outDir)
	if err != nil {
		return nil, err
	}
	if !ok {
		return nil, ErrInvalidIndexDir
	}

	// Mask file
	fileMask := filepath.Join(outDir, MaskFile)
	ok, err = pathutil.Exists(fileMask)
	if err != nil {
		return nil, err
	}
	if !ok {
		return nil, ErrInvalidIndexDir
	}

	// ID list file
	fileIDList := filepath.Join(outDir, IDListFile)
	ok, err = pathutil.Exists(fileIDList)
	if err != nil {
		return nil, err
	}
	if !ok {
		return nil, ErrInvalidIndexDir
	}

	// Trees
	dirTrees := filepath.Join(outDir, TreeDir)
	ok, err = pathutil.DirExists(dirTrees)
	if err != nil {
		return nil, err
	}
	if !ok {
		return nil, ErrInvalidIndexDir
	}

	// ------------- parsing -----------

	idx := &Index{}

	// Mask file
	idx.lh, err = lexichash.NewFromFile(fileMask)
	if err != nil {
		return nil, err
	}
	idx.k = uint8(idx.lh.K)

	// ID list file
	err = idx.readIDlist(fileIDList)
	if err != nil {
		return nil, err
	}

	// Trees
	nMasks := len(idx.lh.Masks)
	idx.Trees = make([]*tree.Tree, nMasks)

	if threads <= 0 {
		threads = runtime.NumCPU()
	}

	treePaths := make([]string, 0, nMasks)
	fs.WalkDir(os.DirFS(dirTrees), ".", func(p string, d fs.DirEntry, err error) error {
		if filepath.Ext(p) == TreeFileExt {
			treePaths = append(treePaths, filepath.Join(dirTrees, p))
		}
		return nil
	})
	if len(treePaths) != nMasks {
		return nil, ErrTreeFileMissing
	}

	var wg sync.WaitGroup
	tokens := make(chan int, threads)
	for _, file := range treePaths {
		wg.Add(1)
		tokens <- 1
		go func(file string) {
			defer func() {
				wg.Done()
				<-tokens
			}()

			// idx of tree
			base := filepath.Base(file)
			i, _err := strconv.Atoi(base[0 : len(base)-len(TreeFileExt)])
			if _err != nil {
				err = _err
				return
			}

			t, _err := tree.NewFromFile(file)
			if _err != nil {
				err = _err
				return
			}

			idx.Trees[i] = t
		}(file)
	}
	wg.Wait()

	if err != nil {
		return nil, err
	}

	for i, t := range idx.Trees {
		if t == nil {
			return nil, fmt.Errorf("tree missing: %d", i)
		}
	}

	return idx, nil
}

func (idx *Index) writeIDlist(file string) error {
	outfh, err := xopen.Wopen(file)
	if err != nil {
		return err
	}
	defer outfh.Close()

	for _, id := range idx.IDs {
		outfh.Write(id)
		outfh.WriteByte('\n')
	}

	return nil
}

func (idx *Index) readIDlist(file string) error {
	fh, err := xopen.Ropen(file)
	if err != nil {
		return err
	}
	defer fh.Close()

	if idx.IDs == nil {
		idx.IDs = make([][]byte, 0, 1024)
	} else {
		idx.IDs = idx.IDs[:0]
	}

	scanner := bufio.NewScanner(fh)
	for scanner.Scan() {
		idx.IDs = append(idx.IDs, []byte(scanner.Text()))
	}
	return scanner.Err()
}
