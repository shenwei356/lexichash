# Changelog

### v0.3.0 - 2024-03-21

- Add 2 methods (`NewWithMasks` and `NewFromTextFile`) for creating LexicHash with custom masks.

### v0.2.0 - 2024-01-29

- Remove all code related to seed indexing and querying and minimize the dependencies.
- Mask() result: only use the last 1 bit, rather than 2, for storing the strand information.
- Add MaskLongSeqs(), which is much faster than mask() for long sequences, like bacteria genomes, requiring nMasks >= 1024.

### v0.1.0 - 2024-01-11

- first fully tested and optimized version.
- Seed indexing and querying are performed in RAM.
