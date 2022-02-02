This repo [implements](lv89.c) the [Landau-Vishkin algorithm][lv89] to compute
the edit distance between two strings. The time complexity of Landau-Vishkin is
O(max(_m_,_n_)+_d_<sup>2</sup>), where _m_ and _n_ are the lengths of the two
sequences and _d_ is their edit distance. This algorithm is faster for more
similar strings.

The actual implementation follows a simplified [WFA][WFA] formulation rather
than the original formulation. It also learns performance tricks from WGA. For
a pair of ~5Mb HLA sequences with ~123k edits, the implementation here can find
the result in 71 seconds, faster than [edlib][edlib].

[lv89]: https://doi.org/10.1016/0196-6774(89)90010-2
[edlib]: https://github.com/Martinsos/edlib
[WFA]: https://github.com/smarco/WFA
