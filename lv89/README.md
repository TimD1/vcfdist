# VCFDist
This project computes the edit distance between two VCF files, when provided
with a reference FASTA.

#### Related Work
This work is based on the [Landau-Vishkin algorithm][lv89] for edit distance,
and follows a simple implementation of [Wave Front Alignment][WFA]. The actual
alignment code is an extension of [Heng Li's implementation][lh3], designed for
computing the edit distance between two highly-similar strings.

[lv89]: https://doi.org/10.1016/0196-6774(89)90010-2
[WFA]: https://github.com/smarco/WFA
[lh3]: https://github.com/lh3/lv89
