# Changelog

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
