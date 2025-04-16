### Summary
Variant correctness within each supercluster is determined by aligning each truth haplotype sequence to a graph constructed from the reference sequence and all query variants. 
By analyzing the path of the optimal alignment (determined by minimum edit distance), we can determine query variant call correctness.
Variants not present on the path are false positives.

We determine which query variants are false positives by seeing if the chosen path through the query graph included or excluded a particular variant. 
We can then assign credit to the remaining variants by marking the points at which all optimal paths diverge and coalesce. We call these "sync points". Between sync points, all query and truth variants are given the same `credit` based on the edit distance of the original truth variants and the current path between the two enclosing sync points. A diagram of this process is shown below. 

<p align="center">
<img src="https://github.com/TimD1/vcfdist/assets/13918078/da1f6d24-ae47-4a45-9be4-9235827a9b4f" alt="precision-recall algorithm"/>
</p>

Some of the minor implementation details of this algorithm may change, but I expect the high-level process to remain the same.


### Derivations for Graph Wavefront Alignment
Below is some basic math for calculating the correct diagonal and offset for the next cell, resulting from a movement into a new adjacent submatrix (for forwards alignment, and backwards path parsing).

#### BACKWARDS QUERY MOVEMENT MATH
```
FORMULA
row = off
col = diag + off
diag = d + 1 - qseqs[qni].length

SITUATION
row = 0

DERIVED
off = 0
col = diag

PREV FORMULA
prev_row = prev_off
prev_col = prev_diag + prev_off
prev_diag = prev_d + 1 - qseqs[prev_qni].length

PREV SITUATION
prev_col = col
prev_row = qseqs[prev_qni].length - 1

PREV DERIVED
prev_col = diag
prev_off = qseqs[prev_qni].length - 1
prev_diag = prev_col - prev_off = diag - qseqs[prev_qni].length + 1
```


#### QUERY MOVEMENT MATH
```
FORMULA
row = off
col = diag + off
diag = d + 1 - qseqs[qni].length

SITUATION
row = qseqs[qni].length - 1

DERIVED
off = qseqs[qni].length - 1
col = d


NEXT FORMULA
next_row = next_off
next_col = next_diag + next_off
next_diag = next_d + 1 - qseqs[next_qni.length]

NEXT SITUATION
next_row = 0
next_col = col

NEXT DERIVED
next_off = 0
next_col = next_diag = d
next_d = d + qseqs[next_qni.length] - 1
```


#### TRUTH MOVEMENT MATH
```
FORMULA
row = off
col = diag + off
diag = d + 1 - qseqs[qni].length

SITUATION
col = tseqs[tni].length - 1


NEXT FORMULA
next_row = next_off
next_col = next_diag + next_off
next_diag = next_d + 1 - qseqs[next_qni.length]

NEXT SITUATION
next_row = row
next_col = 0

NEXT DERIVED
next_off = row
next_diag = next_col - next_off = -row = -off
next_d = next_diag + qseqs[next_qni.length] - 1
next_d = qseqs[next_qni.length] - off - 1

d + 1 - qseqs[qni].length + off = tseqs[tni].length - 1
d + off + 1 = mat_len

diag - d = next_diag - next_d
```


#### BACKWARDS TRUTH MOVEMENT MATH
```
FORMULA
row = off
col = diag + off
diag = d + 1 - qseqs[qni].length

SITUATION
col = 0

DERIVED
diag = -off = -row

PREV FORMULA
prev_row = prev_off
prev_col = prev_diag + prev_off
prev_diag = prev_d + 1 - qseqs[prev_qni.length]

PREV_SITUATION
prev_row = row
prev_col = tseqs[prev_tni].length - 1

PREV DERIVED
prev_off = prev_row = row = -diag
prev_diag = prev_col - prev_off
          = tseqs[prev_tni].length - 1 + diag
prev_d = prev_diag + qseqs[prev_qni].length - 1
```
