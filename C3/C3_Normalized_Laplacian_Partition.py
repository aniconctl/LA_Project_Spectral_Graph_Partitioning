"""
C3: Repeat the Gerrymandering partition with the NORMALIZED Laplacian and compare it with C1.

Same 16-block census graph as C1 (block IDs 2000 through 2015).
    Standard Laplacian    L      = D - A                      (used in C1)
    Normalized Laplacian  L_norm = T L T,  T = D^(-1/2)       (book, Section 1.3)

Three runs are compared. Each run splits the full graph into G1 / G2 and then
splits G1 again into G1a / G1b (the same two-step procedure as Figure 1.23).
    C1            standard Laplacian,   split by the SIGN of the Fiedler vector
    C3            normalized Laplacian, split by the SIGN of the Fiedler vector
    C3 (sweep)    normalized Laplacian, split at the best threshold (smallest Ncut)

The third run is a small extension: it explains why C1 and C3 differ.
Run:  python C3_Normalized_Laplacian_Partition.py   (needs numpy, networkx, matplotlib)
"""

from pathlib import Path
import csv
import numpy as np
import networkx as nx
import matplotlib

matplotlib.use("Agg")  # save figures to files, no window needed
import matplotlib.pyplot as plt

# ---------------------------------------------------------------------------
# 1. The graph (copied from C1: edge list of the 16 census blocks)
# ---------------------------------------------------------------------------
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
A = np.zeros((len(BLOCKS), len(BLOCKS)))
for block_i, block_j in EDGE_PAIRS:
    A[index[block_i], index[block_j]] = A[index[block_j], index[block_i]] = 1

output_dir = Path(__file__).resolve().parent


# ---------------------------------------------------------------------------
# 2. Laplacians and Fiedler vector
# ---------------------------------------------------------------------------
def laplacian(adj):
    """Standard Laplacian L = D - A."""
    return np.diag(adj.sum(axis=1)) - adj


def normalized_laplacian(adj):
    """Normalized Laplacian L_norm = T L T with T = D^(-1/2) (entries 1/sqrt(d_i))."""
    degrees = adj.sum(axis=1)
    if (degrees == 0).any():
        raise ValueError("A vertex has degree 0, so D^(-1/2) does not exist.")
    T = np.diag(1.0 / np.sqrt(degrees))
    return T @ laplacian(adj) @ T


def fiedler_eigenpair(matrix):
    """Return lambda_2 and its unit eigenvector (first listed vertex made positive, as in C1)."""
    eigenvalues, eigenvectors = np.linalg.eigh(matrix)
    vector = eigenvectors[:, 1].copy()
    if vector[0] < 0:
        vector *= -1
    return float(eigenvalues[1]), vector


# ---------------------------------------------------------------------------
# 3. Two ways to cut the Fiedler vector into two groups
# ---------------------------------------------------------------------------
def ncut_value(adj, subset):
    """Ncut of splitting a (sub)graph into `subset` and the rest: cut/vol(S) + cut/vol(rest)."""
    inside = np.zeros(len(adj), dtype=bool)
    inside[subset] = True
    cut = adj[inside][:, ~inside].sum()
    degrees = adj.sum(axis=1)
    return cut / degrees[inside].sum() + cut / degrees[~inside].sum()


def sign_rule(values):
    """Sign rule: values >= 0 form the first group, values < 0 the second group."""
    first = [i for i, v in enumerate(values) if v >= 0]
    second = [i for i, v in enumerate(values) if v < 0]
    return first, second


def sweep_rule(adj, values):
    """Sweep rule: sort the vertices by value, try every threshold, keep the smallest Ncut."""
    order = np.argsort(-values)  # largest value first
    best_k = min(range(1, len(values)), key=lambda k: ncut_value(adj, order[:k]))
    first = sorted(order[:best_k].tolist())
    second = sorted(order[best_k:].tolist())
    return first, second


def split_graph(adj, normalized, rule):
    """Split one (sub)graph in two. Returns lambda_2, Fiedler vector y, scaled vector x, two groups.

    For the normalized Laplacian we also compute x = D^(-1/2) y. It solves L x = lambda D x.
    Since D^(-1/2) has positive entries, x and y always have the same signs.
    """
    degrees = adj.sum(axis=1)
    matrix = normalized_laplacian(adj) if normalized else laplacian(adj)
    lambda2, y = fiedler_eigenpair(matrix)
    x = y / np.sqrt(degrees) if normalized else y.copy()
    first, second = sign_rule(y) if rule == "sign" else sweep_rule(adj, x)
    if 0 not in first:  # name the group that contains the first vertex "first"
        first, second = second, first
    return lambda2, y, x, first, second


def two_step_partition(normalized, rule):
    """Full graph -> G1 / G2, then G1 -> G1a / G1b (induced subgraph, degrees recomputed)."""
    lambda_full, y_full, x_full, first, second = split_graph(A, normalized, rule)
    g1 = [BLOCKS[i] for i in first]
    g2 = [BLOCKS[i] for i in second]

    g1_indices = [index[b] for b in g1]
    g1_adjacency = A[np.ix_(g1_indices, g1_indices)]
    lambda_g1, y_g1, x_g1, first, second = split_graph(g1_adjacency, normalized, rule)
    g1a = [g1[i] for i in first]
    g1b = [g1[i] for i in second]

    return {
        "lambda_full": lambda_full, "y_full": y_full, "x_full": x_full, "G1": g1, "G2": g2,
        "lambda_g1": lambda_g1, "y_g1": y_g1, "x_g1": x_g1, "G1a": g1a, "G1b": g1b,
    }


# ---------------------------------------------------------------------------
# 4. Quality measures of a partition (all computed on the FULL graph)
# ---------------------------------------------------------------------------
DEGREES = A.sum(axis=1)


def partition_metrics(parts):
    """cut edges, RatioCut and Ncut of a list of parts that covers all blocks.

    RatioCut = sum_i cut(S_i, rest) / |S_i|          (divides by the number of blocks)
    Ncut     = sum_i cut(S_i, rest) / vol(S_i)       (divides by the total degree)
    For two parts these equal cut*(1/|S1|+1/|S2|) and cut*(1/vol1+1/vol2).
    """
    ratio_cut = ncut = cut_sum = 0.0
    for part in parts:
        inside = np.array([b in part for b in BLOCKS])
        cut = A[inside][:, ~inside].sum()
        cut_sum += cut
        ratio_cut += cut / inside.sum()
        ncut += cut / DEGREES[inside].sum()
    return {"sizes": [len(p) for p in parts], "cut_edges": int(cut_sum / 2),
            "ratiocut": ratio_cut, "ncut": ncut}


def best_two_way_partitions():
    """Brute force over all 2^15 - 1 ways to cut the 16 blocks into two non-empty groups."""
    n = len(BLOCKS)
    masks = np.arange(2 ** (n - 1) - 1)
    X = np.hstack([np.ones((len(masks), 1)), (masks[:, None] >> np.arange(n - 1)) & 1])
    sizes = X.sum(axis=1)
    volumes = X @ DEGREES
    cuts = np.einsum("ki,ij,kj->k", X, A, 1 - X)
    ratio = cuts * (1 / sizes + 1 / (n - sizes))
    ncut = cuts * (1 / volumes + 1 / (DEGREES.sum() - volumes))
    result = {}
    for name, values in (("best RatioCut", ratio), ("best Ncut", ncut)):
        k = int(np.argmin(values))
        part = [b for b, v in zip(BLOCKS, X[k]) if v == 1]
        rest = [b for b in BLOCKS if b not in part]
        result[name] = (part, rest)
    return result


# ---------------------------------------------------------------------------
# 5. Run everything
# ---------------------------------------------------------------------------
RUNS = {
    "C1": two_step_partition(normalized=False, rule="sign"),
    "C3": two_step_partition(normalized=True, rule="sign"),
    "C3 sweep": two_step_partition(normalized=True, rule="sweep"),
}
best = best_two_way_partitions()

two_way = {name: [run["G1"], run["G2"]] for name, run in RUNS.items()}
two_way["Optimum Ncut (brute force)"] = list(best["best Ncut"])
two_way["Optimum RatioCut (brute force)"] = list(best["best RatioCut"])
three_way = {name: [run["G1a"], run["G1b"], run["G2"]] for name, run in RUNS.items()}

graph = nx.Graph()
graph.add_nodes_from(BLOCKS)
graph.add_edges_from(EDGE_PAIRS)
best_ncut_connected = all(nx.is_connected(graph.subgraph(p)) for p in best["best Ncut"])

# ---------------------------------------------------------------------------
# 6. Save tables
# ---------------------------------------------------------------------------
with (output_dir / "c3_normalized_laplacian.csv").open("w", newline="", encoding="utf-8") as file:
    writer = csv.writer(file)
    writer.writerow(["block_id", *BLOCKS])
    for block, row in zip(BLOCKS, normalized_laplacian(A)):
        writer.writerow([block, *[f"{v:.3f}" for v in row]])

with (output_dir / "c3_fiedler_results.csv").open("w", newline="", encoding="utf-8") as file:
    writer = csv.writer(file)
    writer.writerow(["run", "graph", "block_id", "lambda_2", "fiedler_component", "cluster"])
    for name, run in RUNS.items():
        for block, value in zip(BLOCKS, run["y_full"]):
            cluster = "G1" if block in run["G1"] else "G2"
            writer.writerow([name, "full_graph", block, f"{run['lambda_full']:.3f}",
                             f"{value:.3f}", cluster])
        for block, value in zip(run["G1"], run["y_g1"]):
            cluster = "G1a" if block in run["G1a"] else "G1b"
            writer.writerow([name, "g1_subgraph", block, f"{run['lambda_g1']:.3f}",
                             f"{value:.3f}", cluster])

with (output_dir / "c3_comparison_metrics.csv").open("w", newline="", encoding="utf-8") as file:
    writer = csv.writer(file)
    writer.writerow(["partition", "number_of_clusters", "clusters", "sizes",
                     "cut_edges", "ratiocut", "ncut"])
    for group in (two_way, three_way):
        for name, parts in group.items():
            m = partition_metrics(parts)
            writer.writerow([name, len(parts), " | ".join(" ".join(map(str, p)) for p in parts),
                             "-".join(map(str, m["sizes"])), m["cut_edges"],
                             f"{m['ratiocut']:.3f}", f"{m['ncut']:.3f}"])


# ---------------------------------------------------------------------------
# 7. Figures
# ---------------------------------------------------------------------------
POSITIONS = nx.spring_layout(graph, seed=3, k=0.9, iterations=300)  # same layout in every panel
COLORS = {"G1": "#4C78A8", "G2": "#F58518", "G1a": "#4C78A8", "G1b": "#54A24B"}


def draw_partition(ax, parts, names, title, compare_with=None):
    """Colour the blocks by cluster, draw cut edges dashed red, ring the blocks that moved."""
    label_of = {b: n for p, n in zip(parts, names) for b in p}
    cut = [(i, j) for i, j in EDGE_PAIRS if label_of[i] != label_of[j]]
    kept = [e for e in EDGE_PAIRS if e not in cut]
    nx.draw_networkx_edges(graph, POSITIONS, edgelist=kept, ax=ax, edge_color="#777777", width=1.2)
    nx.draw_networkx_edges(graph, POSITIONS, edgelist=cut, ax=ax, edge_color="#D62728",
                           style="dashed", width=1.8)
    moved = [b for b in BLOCKS if compare_with and label_of[b] != compare_with[b]]
    nx.draw_networkx_nodes(graph, POSITIONS, nodelist=BLOCKS, ax=ax, node_size=620,
                           node_color=[COLORS[label_of[b]] for b in BLOCKS],
                           edgecolors=["black" if b in moved else "white" for b in BLOCKS],
                           linewidths=[3.2 if b in moved else 1.0 for b in BLOCKS])
    nx.draw_networkx_labels(graph, POSITIONS, ax=ax, font_size=7, font_color="white",
                            font_weight="bold")
    m = partition_metrics(parts)
    ax.set_title(f"{title}\n{len(cut)} cut edges, Ncut = {m['ncut']:.3f}", fontsize=10)
    ax.axis("off")
    return label_of


def make_figure(filename, columns, suptitle):
    fig, axes = plt.subplots(1, len(columns), figsize=(5.2 * len(columns), 5.2))
    reference = None
    for ax, (title, parts, names) in zip(axes, columns):
        label_of = draw_partition(ax, parts, names, title, compare_with=reference)
        reference = reference or label_of  # later panels ring blocks that differ from panel 1
    fig.suptitle(suptitle, fontsize=12)
    fig.text(0.5, 0.01, "dashed red = cut edge;  black ring = block placed differently than in C1",
             ha="center", fontsize=9)
    fig.tight_layout(rect=(0, 0.04, 1, 0.95))
    fig.savefig(output_dir / filename, dpi=200)
    plt.close(fig)


r1, r3, rs = RUNS["C1"], RUNS["C3"], RUNS["C3 sweep"]
make_figure("c3_figure_two_clusters.png",
            [("C1: standard Laplacian, sign rule", [r1["G1"], r1["G2"]], ["G1", "G2"]),
             ("C3: normalized Laplacian, sign rule", [r3["G1"], r3["G2"]], ["G1", "G2"]),
             ("C3 sweep: normalized, best threshold", [rs["G1"], rs["G2"]], ["G1", "G2"])],
            "Two clusters of the 16 census blocks")
make_figure("c3_figure_three_clusters.png",
            [("C1: standard Laplacian", [r1["G1a"], r1["G1b"], r1["G2"]], ["G1a", "G1b", "G2"]),
             ("C3: normalized, sign rule", [r3["G1a"], r3["G1b"], r3["G2"]], ["G1a", "G1b", "G2"]),
             ("C3 sweep: normalized, best threshold", [rs["G1a"], rs["G1b"], rs["G2"]],
              ["G1a", "G1b", "G2"])],
            "Three clusters (G1 split again, as in Figure 1.23)")

# Fiedler vectors side by side: which blocks are close to zero?
fig, axes = plt.subplots(2, 1, figsize=(10, 6.2), sharex=True)
for ax, (title, run, key) in zip(axes, [("C1: Fiedler vector of L", r1, "y_full"),
                                        ("C3: Fiedler vector of L_norm", r3, "y_full")]):
    values = run[key]
    ax.bar([str(b) for b in BLOCKS], values,
           color=[COLORS["G1"] if v >= 0 else COLORS["G2"] for v in values])
    ax.axhline(0, color="black", linewidth=0.8)
    ax.set_title(title, fontsize=10)
    ax.set_ylabel("component")
    for b, v in zip(BLOCKS, values):
        if abs(v) < 0.075:
            ax.annotate(f"{v:+.3f}", (str(b), v), ha="center",
                        va="bottom" if v >= 0 else "top", fontsize=8)
axes[1].set_xlabel("census block ID (blue: >= 0, orange: < 0)")
fig.tight_layout()
fig.savefig(output_dir / "c3_figure_fiedler_vectors.png", dpi=200)
plt.close(fig)

# ---------------------------------------------------------------------------
# 8. Print a summary
# ---------------------------------------------------------------------------
for name, run in RUNS.items():
    print(f"\n=== {name} ===")
    print(f"lambda_2 (full graph) = {run['lambda_full']:.3f},  "
          f"lambda_2 (G1 subgraph) = {run['lambda_g1']:.3f}")
    print("G1 :", run["G1"])
    print("G2 :", run["G2"])
    print("G1a:", run["G1a"])
    print("G1b:", run["G1b"])
print("\n=== Two clusters ===")
for name, parts in two_way.items():
    m = partition_metrics(parts)
    print(f"{name:32s} sizes {m['sizes']}  cut {m['cut_edges']}  "
          f"RatioCut {m['ratiocut']:.3f}  Ncut {m['ncut']:.3f}")
print("Best-Ncut split has two connected parts:", best_ncut_connected)
print("\n=== Three clusters ===")
for name, parts in three_way.items():
    m = partition_metrics(parts)
    print(f"{name:32s} sizes {m['sizes']}  cut {m['cut_edges']}  "
          f"RatioCut {m['ratiocut']:.3f}  Ncut {m['ncut']:.3f}")
print("\nWrote tables and figures to", output_dir)
