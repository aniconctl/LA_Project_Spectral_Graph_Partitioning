# C2. Extension: Population Balance Analysis and Improvement

## 1. Assigned Population for the 16 Census Blocks (2000–2015)

The unweighted Fiedler partition in **Figure 1.22** relies purely on graph topology and ignores voter population, violating the legal principle of **equity** (equal population across districts). To evaluate and improve population balance, realistic synthetic population values are assigned to the 16 census blocks (`2000`–`2015`), yielding a total population of $P_{\text{total}} = 2000$ residents and an ideal district population of $P_{\text{ideal}} = 1000$ for $k = 2$ districts.

| Census Block | Baseline District (Fig. 1.22) | Boundary Block? | Assigned Population ($p_i$) |
| :---: | :---: | :---: | :---: |
| 2000 | $G_1$ | No | 130 |
| 2001 | $G_1$ | No | 80 |
| 2002 | $G_1$ | No | 75 |
| 2003 | $G_1$ | No | 115 |
| 2004 | $G_1$ | **Yes** (adj. to 2005, 2006) | 260 |
| 2005 | $G_2$ | **Yes** (adj. to 2004, 2010) | 155 |
| 2006 | $G_2$ | **Yes** (adj. to 2004) | 135 |
| 2007 | $G_2$ | No | 140 |
| 2008 | $G_2$ | No | 160 |
| 2009 | $G_2$ | No | 130 |
| 2010 | $G_1$ | **Yes** (adj. to 2005) | 180 |
| 2011 | $G_1$ | No | 110 |
| 2012 | $G_1$ | No | 55 |
| 2013 | $G_1$ | No | 125 |
| 2014 | $G_1$ | No | 80 |
| 2015 | $G_1$ | No | 70 |
| **Total** | **16 blocks** | — | **2000** |

> **Erratum note:** Page 35 of the textbook lists `2016` inside $G_1$ by mistake. Since the study graph only covers census blocks `2000`–`2015`, block `2016` is omitted here.

---

## 2. Population Balance Measures

For a 2-district partition $(G_1, G_2)$ with ideal population $P_{\text{ideal}} = \frac{P_{\text{total}}}{2} = 1000$:

1. **Maximum Relative Deviation ($\delta_{\max}$):**
   $$\delta_{\max} = \max_{k \in \{1, 2\}} \frac{|P(G_k) - P_{\text{ideal}}|}{P_{\text{ideal}}} \times 100\%$$
2. **Overall Population Range / Total Deviation ($\Delta$):**
   $$\Delta = \frac{\max_k P(G_k) - \min_k P(G_k)}{P_{\text{ideal}}} \times 100\%$$

In U.S. state and local redistricting practice, an overall deviation of $\Delta \le 10\%$ is commonly used as an equity benchmark.

---

## 3. Baseline District Totals (Before Adjustment)

Based on the Fiedler partition in **Figure 1.22**:
- **$G_1$ (11 blocks):** `{2000, 2001, 2002, 2003, 2004, 2010, 2011, 2012, 2013, 2014, 2015}` $\rightarrow P(G_1) = 1280$
- **$G_2$ (5 blocks):** `{2005, 2006, 2007, 2008, 2009}` $\rightarrow P(G_2) = 720$
- **Absolute Difference:** $|1280 - 720| = 560$
- **Balance Metrics:** $\delta_{\max} = 28.000\%$, $\Delta = 56.000\%$ (severely imbalanced; exceeds the $10\%$ benchmark).
- **Cut Size:** $3$ edges — `(2004, 2005)`, `(2004, 2006)`, `(2010, 2005)`.

---

## 4. Improvement Method: Contiguity-Preserving Single Boundary-Block Move

Because $P(G_1) > P(G_2)$, we evaluate moving a single boundary block from the overpopulated district $G_1$ to $G_2$ while verifying that both resulting subgraphs remain connected (**contiguity**). The only boundary blocks in $G_1$ adjacent to $G_2$ are **`2010`** and **`2004`**.

| Partition Scenario | $G_1$ Pop. (Blocks) | $G_2$ Pop. (Blocks) | $|P(G_1)-P(G_2)|$ | $\delta_{\max}$ (%) | $\Delta$ (%) | Cut Edges | Contiguous? |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **Baseline (Fig. 1.22)** | 1280 (11) | 720 (5) | 560 | 28.000 | 56.000 | 3 | Yes |
| Move block `2010` ($G_1 \to G_2$) | 1100 (10) | 900 (6) | 200 | 10.000 | 20.000 | 4 | Yes |
| **Move block `2004` ($G_1 \to G_2$)** | **1020 (10)** | **980 (6)** | **40** | **2.000** | **4.000** | **5** | **Yes** |

### Conclusion
Reassigning boundary block **`2004`** from $G_1$ to $G_2$ produces the most balanced contiguous districts:
- **$G_1'$ (10 blocks):** `{2000, 2001, 2002, 2003, 2010, 2011, 2012, 2013, 2014, 2015}` $\rightarrow P(G_1') = 1020$
- **$G_2'$ (6 blocks):** `{2004, 2005, 2006, 2007, 2008, 2009}` $\rightarrow P(G_2') = 980$

This reduces the overall population deviation $\Delta$ from **56.000%** down to **4.000%** (well within the 10% threshold), at the modest trade-off of increasing the graph cut size from 3 to 5 edges.