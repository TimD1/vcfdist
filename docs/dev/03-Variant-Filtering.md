### Variant Locations
#### BED Region
The `--bed FILENAME` option can be used to select a region to evaluate. Variants outside this region are discarded from the truth and query VCFs.

Variants on the border of BED regions are excluded. This was done to be consistent with Truvari and how most benchmarking truth BEDs were generated (discussion [here](https://github.com/ACEnglish/truvari/issues/193)). In the example below, the BED region is `ref 10 20`. As you can see, Truvari and vcfdist exclude deletions overlapping the border (including if the preceding reference base in the VCF overlaps) whereas vcfeval does not. All three tools exclude insertions on the BED region border.

<p align="center">
<img src="https://github.com/TimD1/vcfdist/assets/13918078/f9d8b9bf-a588-4683-854e-d74dd0ea3f71" alt="bed region example"/>
</p>

#### Spanning deletions
If a variant call is present within a spanning deletion (denoted by an `ALT` field of `*`), it is discarded.

#### Overlapping variants
If two variants overlap or two insertions occur at the exact same position (on the same haplotype), the second variant (and third, etc. if applicable) which would conflict and create an ambiguous or invalid sequence is discarded.

### Variant Attributes

#### Variant Genotype (`GT`)
Variants that are reference calls (`0|0`) or have unknown alleles (`.|.`) are discarded. All other variant genotypes are converted to `1|0` or `0|1` (splitting variants such as `1|1` or `1|2` into two variants and selecting the correct alleles).

#### Genotype separators (`/` vs `|`)
In VCFs, a pipe separator is used for phased variants (`0|1`) and a forward slash is used for unphased variants (`0/1`). All unphased variants are discarded, unless the alleles are the same (`1/1`or `2/2`).

#### Phase set tags (`PS`)
If `PS` is not defined in the VCF header, then all variants are assumed to be globally phased. If a variant is missing a `PS` tag, it is considered to be in the same phase set as adjacent variants missing phase tags (if present).

#### VCF FILTER field
The `--filter FILTER1,FILTER2` option can be used to select variants passing FILTER1 or FILTER2. By default, if this option is not selected, all variants will pass this stage.

#### Variant Size
The `-s SIZE` and `-l SIZE` options can be used to select the smallest and largest variant sizes to be evaluated (inclusive). Other variants are discarded.

#### Variant Quality
The `--min-qual QUAL` and `--max-qual QUAL` options can be used to select the range of variant qualities to evaluate. Variants with lower quality are discarded. Variants with higher quality scores are kept, but their quality score is set to `max_qual`.

### Complex Variants
Complex variants are left and right-trimmed, and then converted to an insertion plus a deletion (or substitution if that is the case).
