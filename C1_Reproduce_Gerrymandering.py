"""
Build the book's 16-block census graph and apply Fiedler splitting
to reproduce the partitions in Figures 1.22 and 1.23.
Rows and columns follow ascending block IDs, 2000 through 2015.
"""

from pathlib import Path
import csv
import numpy as np

"""
Book erratum: the prose for Figure 1.22 includes block 2016,
# but the sample graph contains only block IDs 2000 through 2015.
"""

BLOCKS = list(range(2000, 2016))
EDGE_PAIRS = [
    (2000, 2001), (2000, 2002), (2000, 2003), (2000, 2004),
    (2001, 2002), (2001, 2004), (2001, 2012),
    (2002, 2003), (2002, 2013),
    (2003, 2013),
    (2004, 2005), (2004, 2006), (2004, 2011), (2004, 2012),
    (2005, 2006), (2005, 2009), (2005, 2010),
    (2006, 2007),
    (2007, 2008), (2007, 2009),
    (2010, 2011), (2010, 2014),
    (2011, 2012), (2011, 2013), (2011, 2015),
    (2012, 2013),
    (2013, 2015),
    (2014, 2015),
]

index = {block: i for i, block in enumerate(BLOCKS)}
A = np.zeros((len(BLOCKS), len(BLOCKS)), dtype=int)

for block_i, block_j in EDGE_PAIRS:
    i, j = index[block_i], index[block_j]
    A[i, j] = A[j, i] = 1

degrees = A.sum(axis=1)
L = np.diag(degrees) - A
output_dir = Path(__file__).resolve().parent

def fiedler_eigenpair(laplacian):
    """Return the pair {lambda_2; its unit Fiedler vector}"""
    eigenvalues, eigenvectors = np.linalg.eigh(laplacian)
    vector = eigenvectors[:, 1].copy()

    # Choose the first listed vertex positive.
    if vector[0] < 0:
        vector *= -1

    return float(eigenvalues[1]), vector

def sign_partition(labels, vector, positive_name, negative_name):
    """Split vertices at zero and return groups plus label assignments"""
    positive = [label for label, value in zip(labels, vector) if value >= 0]
    negative = [label for label, value in zip(labels, vector) if value < 0]

    assignments = {label: positive_name for label in positive}
    assignments.update({label: negative_name for label in negative})

    return positive, negative, assignments

# First split: calculate G1 and G2 from the full graph's Fiedler vector
lambda2_full, fiedler_full = fiedler_eigenpair(L)
g1_blocks, g2_blocks, full_assignments = sign_partition(BLOCKS, fiedler_full, "G1", "G2")

# Second-stage split: build the induced subgraph from the calculated G1
g1_indices = [index[block] for block in g1_blocks]
g1_adjacency = A[np.ix_(g1_indices, g1_indices)]

g1_set = set(g1_blocks)
g1_edge_pairs = [(i, j) for i, j in EDGE_PAIRS if i in g1_set and j in g1_set]

# Recompute degrees and the Laplacian within G1
g1_degrees = g1_adjacency.sum(axis=1)
g1_laplacian = np.diag(g1_degrees) - g1_adjacency

lambda2_g1, fiedler_g1 = fiedler_eigenpair(g1_laplacian)
g1a_blocks, g1b_blocks, g1_assignments = sign_partition(g1_blocks, fiedler_g1, "G1a", "G1b")

# Figure 1.23 keeps G2 and splits G1 into G1a and G1b.
figure_1_23_partitions = {"G1a": g1a_blocks, "G1b": g1b_blocks, "G2": g2_blocks}

# Save the full graph and the G1 subgraph matrices
def write_labeled_matrix(filename, labels, matrix):
    """Write a matrix CSV with block IDs as row and column labels"""
    with (output_dir / filename).open("w", newline="", encoding="utf-8") as file:
        writer = csv.writer(file)
        writer.writerow(["block_id", *labels])
        for block_id, row in zip(labels, matrix):
            writer.writerow([block_id, *row.tolist()])

write_labeled_matrix("binary_adjacency_A.csv", BLOCKS, A)
write_labeled_matrix("laplacian_L.csv", BLOCKS, L)
write_labeled_matrix("figure_1_23_g1_subgraph_adjacency.csv", g1_blocks, g1_adjacency)
write_labeled_matrix("figure_1_23_g1_subgraph_laplacian.csv", g1_blocks, g1_laplacian)

with (output_dir / "edge_list.csv").open("w", newline="", encoding="utf-8") as file:
    writer = csv.writer(file)
    writer.writerow(["block_i", "block_j"])
    writer.writerows(EDGE_PAIRS)

with (output_dir / "figure_1_23_g1_subgraph_edges.csv").open("w", newline="", encoding="utf-8") as file:
    writer = csv.writer(file)
    writer.writerow(["block_i", "block_j"])
    writer.writerows(g1_edge_pairs)

with (output_dir / "fiedler_results.csv").open("w", newline="", encoding="utf-8") as file:
    writer = csv.writer(file)
    writer.writerow(["graph", "block_id", "lambda_2", "fiedler_component", "cluster"])

    for block, value in zip(BLOCKS, fiedler_full):
        writer.writerow(
            ["full_graph", block, f"{lambda2_full:.3f}", f"{value:.3f}", full_assignments[block]]
        )

    for block, value in zip(g1_blocks, fiedler_g1):
        writer.writerow(
            ["g1_subgraph", block, f"{lambda2_g1:.3f}", f"{value:.3f}", g1_assignments[block]]
        )

print(f"Full sample: {len(BLOCKS)} vertices, {len(EDGE_PAIRS)} edges")
print("Degree sequence:", degrees.tolist())
print(f"Figure 1.22 Fiedler value (lambda_2): {lambda2_full:.3f}")
print("Figure 1.22 G1:", g1_blocks)
print("Figure 1.22 G2:", g2_blocks)
print("Figure 1.22 Fiedler vector:", dict(zip(BLOCKS, np.round(fiedler_full, 3))))
print(f"Figure 1.23 G1 subgraph Fiedler value (lambda_2): {lambda2_g1:.3f}")
print("Figure 1.23 G1a:", g1a_blocks)
print("Figure 1.23 G1b:", g1b_blocks)
print("Figure 1.23 retained G2:", g2_blocks)
print("Figure 1.23 three-cluster partition:", figure_1_23_partitions)
print(f"G1 subgraph: {len(g1_blocks)} vertices, " f"{len(g1_edge_pairs)} induced edges")
print("Wrote matrices, edge lists, and Fiedler results to", output_dir)