# Randomized-Layered-ADMM-LP
This repository includes the codes for the novel randomized scheduling strategy introduced in the paper:
[Randomized Scheduling of ADMM-LP Decoding Based on Geometric Priors](https://ieeexplore.ieee.org/stamp/stamp.jsp?arnumber=9965857), 
presented in ITW '22, 
and part of my masters thesis,
[Approximate and Randomized ADMM-LP Decoding Using the Geometric Information of the Parity Polytope](https://tspace.library.utoronto.ca/bitstream/1807/125679/1/Asadzadeh_Amirreza_202211_MAS_thesis.pdf).

## Background
This project is focused on decoding low-density parity-check (LDPC) codes using an alternating direction method of multipliers (ADMM) framework
for solving the linear-programming (LP) decoding problem. This algorithm was initially introduced in the paper [Decomposition Methods for Large Scale LP Decoding](https://ieeexplore.ieee.org/abstract/document/6595057) and was significantly improved in the papers
[the ADMM Penalized Decoder for LDPC Codes](https://ieeexplore.ieee.org/abstract/document/7456284) and [Hardware Based Projection onto the Parity Polytope and Probability Simplex](https://ieeexplore.ieee.org/abstract/document/7421292).

ADMM-LP decoding algorithm iteratively applies message passing decoding on the Tanner graph of LDPC codes, while stroing the residual information in the Lagrange multipliers.
As part of this algorithm, in each iteration, there exists $M$ projections onto a geometric object known as parity polytope, for all $M$ check nodes of the graph.
Such projections can be implemented in practice via a water-filling process, which includes sorting and thresholding operations.
The average time complexity of projecting a vector of size $d$ onto its parity polytope of dimension $d$ is $O(d \log(d))$, due to the sorting operation.
It is known that this projection step is the complexity bottleneck of ADMM-LP decoding, and simplifying this step results in the complexity reduction of the overall decoding algorithm.

## Idea
In this project, we introduce a novel randomized scheduling strategy in order to boost the convergence of the iterative ADMM-LP decoding algorithm.
Assume we have $M$ check nodes in the Tanner graph of the LDPC code, each for one parity-check constraint.
We induce a finite probability distribution over the indices of those check nodes, then we choose one check node in each decoding round based on the introduced probability.
Specifically, we induce the [Boltzmann distribution](https://en.wikipedia.org/wiki/Boltzmann_distribution) (or equivalently,
apply the [Softmax function](https://en.wikipedia.org/wiki/Softmax_function)) over the check indices,
such that the probability of choosing the $j$-th check node with replica vector $z_j$ of dimension $d_j$ becomes
$$\pi_j = \frac{\exp(-\phi||z_j - 0.5\mathbf{1}_{d_j}||_2)}{\mathcal{Z}},$$
where $\phi$ and $\mathcal{Z}$ are the hyperparameter of the model and the normalization parameter, respectively.
