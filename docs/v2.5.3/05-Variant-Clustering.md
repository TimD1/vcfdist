In order to make whole-genome alignment-based evaluation of variant calling tractable, we need a method of breaking contig-sized alignments (~100Mb) into many smaller sub-problems (each <50kb). This is done by identifying dependent variant calls through clustering, a process that occurs in two stages. Firstly, variants are clustered within their haplotype (there are 4 in total, from diploid truth and query VCFs). Then, a second stage which we call "superclustering" groups variants across all four haplotypes.

### Clustering
#### Summary
There are currently three implemented clustering methods: biwfa (default, most accurate), gap N (simple, efficient), and size N (efficient, good for large SVs). For each method, the left and right "reaches" are first calculated, which define the genomic interval for which this variant is not independent of other variants. Any variants with overlapping intervals are merged into a single cluster. Only `biwfa` occurs iteratively, with up to `-i` iterations.

#### `biwfa`
```
left_reach = variant_end - [maximum reach of all leftward alignments starting from variant_end 
                            with lesser alignment penalty to current variant(s) representation]
right_reach = variant_start - [maximum reach of all rightward alignments starting from variant_start 
                               with lesser alignment penalty to current variant(s) representation]
```
WFA is a time and space efficient affine-gap alignment algorithm, which we use bi-directionally to find possible alternate variant alignments (and therefore if nearby variants are independent). See the papers on [BiWFA](https://academic.oup.com/bioinformatics/article/39/2/btad074/7030690) and [WFA](https://academic.oup.com/bioinformatics/article/37/4/456/5904262) for more details. This is the currently recommended (and default) vcfdist clustering algorithm because it is the most accurate; it will always find dependencies if they exist. However, when evaluating large structural variants (above 1kbp) it tends to create large clusters (10-50kbp), which results in large memory usage and slower evaluations. For evaluating large variants, using `--cluster size 100` may be preferable.

<p align="center">
<img src="https://github.com/TimD1/vcfdist/assets/13918078/8f2f5b35-5a23-4079-ab21-d1660dc88e6a" alt="biwfa clustering"/>
</p>

#### `gap GAP_WIDTH`
```
left_reach = variant_start - GAP_WIDTH
right_reach = variant_end + GAP_WIDTH
```
Gap-based clustering is the simplest and fastest clustering method: group together all variants less than N bases apart. It is also the least accurate, and will miss variant dependencies if GAP_WIDTH is too small. Conversely, as GAP_WIDTH nears the reciprocal of the background rate of genomic variation between humans (one SNP every 1000 bases), clusters will be merged on average and will grow to be very large. We recommend 50 < GAP_WIDTH < 200, and to limit evaluations to small variants when using this option.

#### `size GAP_WIDTH`
```
left_reach = variant_start - std::max(GAP_WIDTH, sizeof(variant))
right_reach = variant_end + std::max(GAP_WIDTH, sizeof(variant))
```
This is a heuristic that compromises in terms of efficiency and accuracy, basically extending the gap heuristic to work with larger variants. Once a variant is larger than size GAP_WIDTH, the required gap to consider it independent of an adjacent variant is the size of the variant, instead of GAP_WIDTH.

### Superclustering
Using the calculated left and right reaches for clusters on each of the four haplotypes, a merge occurs across all four haplotypes whenever cluster reaches overlap, into a supercluster. If there is no overlap, a cluster is converted into a supercluster nevertheless.

<p align="center">
<img width="800px" src="https://github.com/TimD1/vcfdist/assets/13918078/e1690116-0690-4776-86c9-a5d57d5b95a6" alt="superclustering"/>
</p>
In the above image, clusters within the first supercluster are assumed to have overlapping left and right reaches, even if the variants within the clusters themselves don't always overlap.