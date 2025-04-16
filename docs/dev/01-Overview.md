vcfdist evaluates the correctness of a set of phased variant calls (query VCF) relative to a set of phased ground truth variant calls (truth VCF) for a subset (regions BED) of the desired genome (reference FASTA). vcfdist was designed to evaluate human genomes, but should work on other monoploid and diploid species. It can evaluate variants of any type, including STRs (simple tandem repeats) and CNVs (copy number variants), but vcfdist classifies variants into SNPs (single nucleotide polymorphisms), INDELS (insertions and deletions), and SVs (structural variants) during evaluation. Evaluating variants larger than 10,000 bases is not recommended at the moment, as it will require large amounts of memory (over 50GB RAM). Below is a diagrammatic overview of vcfdist. Inputs are shown in red, internal steps in yellow, and optional steps in gray.

<p align="center"><img src="https://github.com/TimD1/vcfdist/assets/13918078/85ecfbdb-0028-4bc8-aae0-3422849d1fd2" alt="overview"/></p>

## Index
 - [Parameters and Usage](https://github.com/TimD1/vcfdist/wiki/02-Parameters-and-Usage)
 - [Variant Filtering](https://github.com/TimD1/vcfdist/wiki/03-Variant-Filtering)
 - [VCF Normalization](https://github.com/TimD1/vcfdist/wiki/04-VCF-Normalization)
 - [Variant Clustering](https://github.com/TimD1/vcfdist/wiki/05-Variant-Clustering)
 - [Precision and Recall](https://github.com/TimD1/vcfdist/wiki/06-Precision-and-Recall)
 - [Phasing Analysis](https://github.com/TimD1/vcfdist/wiki/07-Phasing-Analysis)
 - [Alignment Distance](https://github.com/TimD1/vcfdist/wiki/08-Alignment-Distance)
 - [Outputs](https://github.com/TimD1/vcfdist/wiki/09-Outputs)
 - [Variant Stratification](https://github.com/TimD1/vcfdist/wiki/10-Variant-Stratification)

## Repository Structure
<table>
<tr>
  <th> Folder </th>
  <th> Description </th>
</tr>
<tr>
  <td> <code>src</code> </td>
  <td> C++ source code for vcfdist</td>
</tr>
<tr>
  <td> <code>demo</code> </td>
  <td> a simple self-contained vcfdist example script, including inputs and expected output</td>
</tr>
<tr>
  <td> <code>analysis</code> </td>
  <td> analysis scripts for "vcfdist: accurately benchmarking phased small variant calls" </td>
</tr>
<tr>
  <td> <code>analysis-v2</code> </td>
  <td> analysis scripts for "Jointly benchmarking phased small and structural variant calls with vcfdist" </td>
</tr>
<tr>
  <td> <code>assets</code> </td>
  <td> project assets for Github and Wiki (images, LaTeX source, etc) </td>
</tr>
<tr>
  <td> <code>docs</code> </td>
  <td> old Wiki documentation </td>
</tr>
</table>
