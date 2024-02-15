#### Phasing a single supercluster
In order to phase a single supercluster, we first perform four separate alignments and store the edit distance (ED) for each.
1. Query Haplotype #1 to Truth Haplotype #1 (`Q1T1`)
2. Query Haplotype #2 to Truth Haplotype #2 (`Q2T2`)
3. Query Haplotype #1 to Truth Haplotype #2 (`Q1T2`)
4. Query Haplotype #2 to Truth Haplotype #1 (`Q2T1`)

Then, the two possible phasings are compared by summing the edit distances for a corresponding pair of alignments.

Phasing A (original): `ED(Q1T1) + ED(Q2T2) = A`

Phasing B  (flipped): `ED(Q1T2) + ED(Q2T1) = B`

If there is a significant difference between the edit distance resulting from the two possible phasings, the supercluster is considered phased:

```
(max(A,B) - min(A,B)) / max(A,B) > phasing_threshold ? phased : unphased;
if (phased) { phasing = A < B ? 'A' : 'B'; }
```

The default threshold is `--phasing-threshold 0.6`. A value of `0.0` will simply choose the phasing with lower edit distance, whereas a value of `1.0` will require one phasing to have an edit distance of zero in order to consider the supercluster phased. A threshold of `0.6` requires at least a 60% reduction in edit distance to consider the supercluster confidently phased.

#### Calculating switch and flip errors
Once per-supercluster phasings are assigned, a dynamic programming algorithm is used to minimize the total number of switch and flip errors within each phase block. Black arrows represent penalized transitions, whereas colored arrows are not penalized. Note that changing the phasing threshold will affect which superclusters are considered phased, thus affecting the number of reported flip errors (an unphased supercluster will never be a flip error).

<p align="center">
<img src="https://github.com/TimD1/vcfdist/assets/13918078/f3d0a11b-322a-47d6-b670-ae854fdab089" alt="phasing algorithm overview"/>
</p>

Note: The label "Supercluster Phasing Error" in the above diagram we now refer to as a supercluster flip error.

#### Metric Definitions
Phase set NG50 is the largest phased region such that all phased regions of that size and larger cover at least 50% of the genome.
Phase block switch NGC50 is calculated similarly, but instead of only breaking regions at new phase sets, it also breaks regions when switch errors occur.
Phase block switchflip NGC50 additionally breaks regions on flip errors.