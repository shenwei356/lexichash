# lexichash

[![Go Reference](https://pkg.go.dev/badge/github.com/shenwei356/lexichash.svg)](https://pkg.go.dev/github.com/shenwei356/lexichash)

This project implements [LexicHash](https://academic.oup.com/bioinformatics/article/39/11/btad652/7329717) in Golang,
with high performance and a low memory footprint.

This package uses custom radix trees to store each mask's bit-packed k-mers with location information,
and indexing, serialization and searching methods are also provided.

## Related projects

- This package is used in [LexicMap](https://github.com/shenwei356/LexicMap).
- Bit-packed k-mer operations is provided by [kmers](https://github.com/shenwei356/kmers).
- The [radix tree data structure](https://github.com/shenwei356/lexichash/tree/main/tree) was originally modified from [go-radix](https://github.com/armon/go-radix),
  then the structure was simplified for k-mers and more methods were added.

## Support

Please [open an issue](https://github.com/shenwei356/lexichash/issues) to report bugs,
propose new functions or ask for help.

## License

[MIT License](https://github.com/shenwei356/lexichash/blob/master/LICENSE)

