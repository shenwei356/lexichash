# Changelog

### v0.6.0 - 2026-09-07

- Change parameter `skipRegions` from `[][2]int` to a flat `[]int` of start/end pairs in all masking methods.
- Add a function `SortSkipRegions` to sort the regions in place.
- Fix `Mask` and `MaskLongSeqs` not applying `skipRegions` in the second round.

### v0.5.5 - 2026-08-21

- Add `MaskKnownDistinctPrefixesWithStrandBias` for compatibility with LexicMap indexes v3.4 and earlier.

### v0.5.4 - 2026-08-20

- **Fix `MaskKnownDistinctPrefixes` incorrectly skipping some negative-strand k-mers.**
- Improve the speed of `MaskKnownPrefixes` and `MaskKnownDistinctPrefixes` by 20%.

### v0.5.3 - 2026-06-02

- Slightly speedup

### v0.5.2 - 2026-03-04

- Check the validity of the sequence; never trust the user.

### v0.5.1 - 2026-02-28

- Added a new method `SupportSoftMasking` for `LexicHash`.

### v0.5.0 - 2024-07-22

- Better mask generation: requiring prefixes of `p+1` to be distinct.
- New methods `IndexMasksWithDistinctPrefixes` and `MaskKnownDistinctPrefixes`.
- Safer `MaskKmer`.

### v0.4.2 - 2024-06-16

- Fix skipping regions by further including the last k-1 bases of contigs.

### v0.4.1 - 2024-06-12

- Faster `MaskKnownPrefixes` by replacing `map` with `slice` to act as the lookup table of prefixes.

### v0.4.0 - 2024-05-13

- All A's or N's are skipped in k-mer generation step.
- Add method `IndexMasks` and `MaskKnownPrefixes` for faster masking k-mers of which the prefixes are existed.
- Add method `MaskKmer` to returns the indexes of masks that possibly mask a k-mer. 

### v0.3.0 - 2024-03-21

- Add 2 methods (`NewWithMasks` and `NewFromTextFile`) for creating LexicHash with custom masks.

### v0.2.0 - 2024-01-29

- Remove all code related to seed indexing and querying and minimize the dependencies.
- Mask() result: only use the last 1 bit, rather than 2, for storing the strand information.
- Add MaskLongSeqs(), which is much faster than mask() for long sequences, like bacteria genomes, requiring nMasks >= 1024.

### v0.1.0 - 2024-01-11

- first fully tested and optimized version.
- Seed indexing and querying are performed in RAM.
