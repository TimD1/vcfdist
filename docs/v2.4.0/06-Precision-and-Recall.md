### Summary
Variant correctness is determined by first calculating the minimum edit distance alignment of the query "sequence" to the truth sequence. The query "sequence" is not quite a sequence but instead a specific type of graph which allows arbitrary omission of query variants (so that we can ignore false positives). By backtracking through all optimal alignments, we can determine which query variants are false positives by seeing if the chosen path through the query graph included or excluded a particular variant. We can then assign credit to the remaining variants by marking the points at which all optimal paths diverge and coalesce. We call these "sync points". Between sync points, all query and truth variants are given the same `credit` based on the edit distance of the original truth variants and the current path between the two enclosing sync points. A diagram of this process is shown below. 

<p align="center">
<img src="https://github.com/TimD1/vcfdist/assets/13918078/da1f6d24-ae47-4a45-9be4-9235827a9b4f" alt="precision-recall algorithm"/>
</p>

Some of the minor implementation details of this algorithm may change, but I expect the high-level process to remain the same.

