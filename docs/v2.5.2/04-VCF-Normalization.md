Following variant [clustering](https://github.com/TimD1/vcfdist/wiki/05-Variant-Clustering), variants are optionally realigned by selecting the `--realign-query` and/or `--realign-truth` options. The `--realign-only` flag can be used to skip downstream evaluations.

### Best-Alignment Normalization
As initially introduced in [this manuscript](https://doi.org/10.1093/bioinformatics/btw748) and further explored in [our work](https://doi.org/10.1038/s41467-023-43876-x), best alignment normalization can be used to select between several possible variant representations when complex variants are involved. Affine gap [Smith Waterman](https://en.wikipedia.org/wiki/Smith%E2%80%93Waterman_algorithm) alignment is used to select the "best" variant representation, defined by a given set of alignment parameters. The design space for these parameters (m, x, o, e) is shown below, with many common alignment tools plotted and four example (A,B,C,D) alignments with their resulting variant  representations. By default, the representation selected by vcfdist is at Point C.
<p align="center">
<img src="https://github.com/TimD1/vcfdist/assets/13918078/047cdde9-57d1-4625-993d-49071e2b6095" alt="the affine-gap design space for variant representation"/>
</p>

### Standard VCF Normalization
The traditional method of variant normalization involves decomposing complex variants, trimming unnecessary bases from the variant representation, and then left-aligning INDELs.
<p align="center">
<img src="https://github.com/TimD1/vcfdist/assets/13918078/bf15f99c-1c00-4abe-bc03-d8d2afff1cf0" alt="variant decomposition, trimming, and left-shifting"/>
</p>
This procedure is sufficient to create a unique canonical representation for a single variant, but not when multiple or complex variants are involved.