import networkx as nx

# 1. Define the 16-block Charlotte census graph (blocks 2000-2015) from Figure 1.22
edges = [
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

G = nx.Graph()
G.add_edges_from(edges)

# 2. Assigned synthetic populations (Total = 2000, Ideal = 1000 per district)
population = {
    2000: 130, 2001: 80,  2002: 75,  2003: 115,
    2004: 260, 2005: 155, 2006: 135, 2007: 140,
    2008: 160, 2009: 130, 2010: 180, 2011: 110,
    2012: 55,  2013: 125, 2014: 80,  2015: 70,
}

# 3. Baseline Fiedler partition (Figure 1.22)
G1 = {2000, 2001, 2002, 2003, 2004, 2010, 2011, 2012, 2013, 2014, 2015}
G2 = {2005, 2006, 2007, 2008, 2009}

def evaluate_partition(name, d1, d2):
    p1 = sum(population[u] for u in d1)
    p2 = sum(population[u] for u in d2)
    p_ideal = (p1 + p2) / 2.0
    delta_max = max(abs(p1 - p_ideal), abs(p2 - p_ideal)) / p_ideal * 100
    delta_range = abs(p1 - p2) / p_ideal * 100
    cut_edges = nx.cut_size(G, d1, d2)
    contiguous = nx.is_connected(G.subgraph(d1)) and nx.is_connected(G.subgraph(d2))
    print(f"{name:32s} | P(G1)={p1:4d} | P(G2)={p2:4d} | "
          f"d_max={delta_max:6.3f}% | Range={delta_range:6.3f}% | "
          f"Cut={cut_edges} | Contiguous={contiguous}")

# Evaluate baseline
evaluate_partition("Baseline (Figure 1.22)", G1, G2)

# Find and test all boundary blocks in the overpopulated district G1
boundary_G1 = sorted(u for u in G1 if any(v in G2 for v in G.neighbors(u)))
for block in boundary_G1:
    evaluate_partition(f"Move block {block} (G1 -> G2)", G1 - {block}, G2 | {block})