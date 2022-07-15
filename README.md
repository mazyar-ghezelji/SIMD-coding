# SIMD-coding
SIMD (single instruction multiple output) implementation of local, global, and semi-global sequence alignment algorithms.

# ABSTRACT
Sequence alignment is one of the most crucial steps in bioinformatics analysis. Many algorithms have been developed with the goal of generating optimal solutions to this problem in an efficient manner.<br />

But, optimality in this context, which is more important in some use cases, comes with higher orders of complexity, meaning that the execution time for large datasets would be a bottleneck. One of the ways to reduce the execution time for these algorithms is to utilize parallelization. <br />

SIMD (single input multiple outputs) instruction-level parallelization is one of the powerful tools that can reduce the run time of the algorithms substantially. There are tools like SEQ programming language that use this tool to increase their efficiency. <br />

In this report I am implementing global, local, and semi-global sequence alignment algorithms using SIMD and we compare the results with SEQ alignment modules.<br />
