# C3: Partition with the normalized Laplacian and comparison with C1

Group project, Topic 1.2 (Spectral Graph Theory and Partitioning), Part C3:
*repeat the Gerrymandering partition with the normalized Laplacian and compare it with C1.*

Everything here is produced by `C3_Normalized_Laplacian_Partition.py`. The graph is the same 16-block census graph as in C1 (blocks 2000 to 2015, 28 edges).

## 1. Files

| File | Content |
|---|---|
| `C3_Normalized_Laplacian_Partition.py` | code: C1 and C3 partitions, quality measures, figures |
| `c3_figure_two_clusters.png` | two clusters: C1 vs C3 vs C3 (sweep) |
| `c3_figure_three_clusters.png` | three clusters (G1 split again, as in Figure 1.23) |
| `c3_figure_fiedler_vectors.png` | Fiedler vector of L and of the normalized Laplacian |
| `c3_fiedler_results.csv` | Fiedler vectors and cluster labels of all three runs |
| `c3_comparison_metrics.csv` | cut edges, RatioCut, Ncut of every partition |
| `c3_normalized_laplacian.csv` | the 16 x 16 normalized Laplacian (3 decimals) |

## 2. What changes with the normalized Laplacian

* C1 uses the standard Laplacian `L = D - A`.
* C3 uses the normalized Laplacian `L_norm = T L T`, where `T` is diagonal with `T_ii = 1/sqrt(d_i)` (book, Section 1.3). Its entries are `1` on the diagonal and `-1/sqrt(d_i d_j)` if blocks `i` and `j` share an edge. Example: blocks 2000 and 2003 have degrees 4 and 3, so the entry is `-1/sqrt(12) = -0.289`.
* The Fiedler vector `y` is the unit eigenvector of the 2nd smallest eigenvalue, and we split by the sign of its entries, exactly as in C1. Then G1 is split again using the induced subgraph of G1, with its degrees recomputed.
* Two facts that explain the difference:
  1. For `L`, the Fiedler vector satisfies `sum_i y_i = 0` (it is orthogonal to the all-ones vector). For `L_norm`, it satisfies `sum_i sqrt(d_i) y_i = 0`. So in C3 each block counts proportionally to its **degree**, and "balanced around zero" now means balanced in total degree.
  2. If `x = D^(-1/2) y`, then `L x = lambda D x`. Because `D^(-1/2)` has positive entries, `x` and `y` have the same signs, so using `x` instead of `y` would not change the sign partition.
* The eigenvalues of `L_norm` are between 0 and 2, so the values of lambda_2 of the two Laplacians cannot be compared with each other. For the full graph: lambda_2(L) = 0.338, lambda_2(L_norm) = 0.146.

## 3. Results

### 3.1 Two clusters (Figure `c3_figure_two_clusters.png`)

| | G1 | G2 |
|---|---|---|
| C1 (standard L) | {2000, 2001, 2002, 2003, 2004, 2010, 2011, 2012, 2013, 2014, 2015} | {2005, 2006, 2007, 2008, 2009} |
| C3 (normalized) | {2000, 2001, 2002, 2003, 2011, 2012, 2013, 2014, 2015} | {2004, 2005, 2006, 2007, 2008, 2009, 2010} |

C3 moves blocks **2004** and **2010** from G1 to G2. Cut edges:

* C1 cuts 3 edges: (2004,2005), (2004,2006), (2005,2010).
* C3 cuts 6 edges: (2000,2004), (2001,2004), (2004,2011), (2004,2012), (2010,2011), (2010,2014).

### 3.2 Three clusters (Figure `c3_figure_three_clusters.png`)

| | G1a | G1b | G2 |
|---|---|---|---|
| C1 | {2000, 2001, 2002, 2003, 2004, 2012, 2013} | {2010, 2011, 2014, 2015} | {2005, 2006, 2007, 2008, 2009} |
| C3 | {2000, 2001, 2002, 2003} | {2011, 2012, 2013, 2014, 2015} | {2004, 2005, 2006, 2007, 2008, 2009, 2010} |

(The C1 row is the partition of the book's Figure 1.23. For the second step, lambda_2 is 0.675 in C1 and 0.277 in C3.)

## 4. Comparison

We measure each partition on the **full graph** with three numbers (smaller is better):

* **cut edges**: number of edges between different clusters;
* **RatioCut** = sum over clusters of `cut(S) / |S|` (divides by the number of blocks);
* **Ncut** = sum over clusters of `cut(S) / vol(S)`, where `vol(S)` is the sum of the degrees in `S`.

Two clusters:

| Partition | Sizes | Cut edges | RatioCut | Ncut |
|---|---|---|---|---|
| C1 (standard L, sign rule) | 11-5 | 3 | 0.873 | 0.301 |
| C3 (normalized, sign rule) | 9-7 | 6 | 1.524 | 0.449 |
| C3 sweep (normalized, best threshold) | 11-5 | 3 | 0.873 | 0.301 |
| Best possible Ncut (all 32,767 splits) | 11-5 | 3 | 0.873 | 0.301 |
| Best possible RatioCut (all 32,767 splits) | 13-3 | 2 | 0.821 | 0.373 |

Three clusters:

| Partition | Sizes | Cut edges | RatioCut | Ncut |
|---|---|---|---|---|
| C1 | 7-4-5 | 7 | 2.707 | 0.815 |
| C3 (normalized, sign rule) | 4-5-7 | 9 | 3.507 | 0.974 |
| C3 sweep (normalized, best threshold) | 7-4-5 | 7 | 2.707 | 0.815 |

### What we observe

1. **C1 and C3 are different partitions.** With two clusters they disagree on 2 blocks out of 16 (2004 and 2010). With three clusters, 4 blocks (2004, 2010, 2012 and 2013) are grouped differently.
2. **In this example C3 (sign rule) is worse than C1** by every measure: 6 cut edges instead of 3 for two clusters, and 9 instead of 7 for three clusters. The cluster sizes are the same set (4, 5 and 7 blocks) in both three-cluster partitions, so balance does not separate them.
3. **The C1 two-cluster partition is the best possible for Ncut.** We checked all 32,767 ways to cut the 16 blocks into two groups: the smallest Ncut is 0.301, reached by the C1 partition (and both its parts are connected). The smallest RatioCut (0.821) is reached by a different split that cuts only 2 edges but leaves only 3 blocks on one side.
4. **Why does C3 miss it?** The sign rule uses the threshold 0, which is not guaranteed to be the best threshold. Two blocks sit almost at zero: in the normalized Fiedler vector block 2004 has value -0.001 and block 2010 has value -0.030 (in C1 they are 0.059 and 0.071). They are boundary blocks, and a tiny change of the threshold moves them. See `c3_figure_fiedler_vectors.png`.
5. **The sweep rule fixes this.** If we keep the normalized Laplacian but, instead of the threshold 0, try every threshold along the sorted entries of `x = D^(-1/2) y` and keep the one with the smallest Ncut, we get **exactly the C1 partition**, for two clusters and for three clusters.

## 5. Conclusion

* The normalized Laplacian gives a different partition when the sign of the Fiedler vector is used with threshold 0, because it weights each block by its degree, and boundary blocks such as 2004 (6 neighbours) and 2010 land on the other side.
* For this small census-block graph, the standard Laplacian (C1) is the better choice with the sign rule. The normalized Laplacian with a best-threshold (sweep) search reaches the same answer as C1.
* Limits: one small unweighted graph with 16 blocks, so this does not show that one Laplacian is better in general. Also, none of these cuts looks at population or shape (equity and compactness), only at adjacency.

## 6. How to run

```
python C3_Normalized_Laplacian_Partition.py
```

Needs `numpy`, `networkx`, `matplotlib`. The files above are written in the same folder as the script. Running it takes a few seconds.

## References

* C. Arangala (2023), *Linear Algebra With Machine Learning and Data*, CRC Press: normalized Laplacian, p. 14; Fiedler vector and the case study, pp. 21 and 34 to 36.
* J. Shi and J. Malik (2000), *Normalized Cuts and Image Segmentation*, IEEE Transactions on Pattern Analysis and Machine Intelligence 22(8), 888-905, https://publications.ri.cmu.edu/normalized-cuts-and-image-segmentation (accessed 22 November 2026).
* U. von Luxburg (2007), *A Tutorial on Spectral Clustering*, Statistics and Computing 17(4), https://arxiv.org/abs/0711.0189 (accessed 22 November 2026).
