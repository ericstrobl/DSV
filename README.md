# Differentially Supervised Varimax (DSV)

Most component-analytic methods (e.g., PCA, factor analysis, standard varimax) prioritize shared variance among behavioral items and therefore recover dimensions that summarize the overall covariance structure. They do not, by design, isolate between-group differences in how behaviors relate to clinically meaningful distress dimensions.

Differentially Supervised Varimax (DSV) addresses a complementary inferential target: it asks which item-defined behavioral profiles show systematically different associations with distress dimensions across two groups (e.g., ASD vs TD), while preserving orthogonality and interpretability.

DSV takes: 

- `X`: specific / clinically interpretable severity dimensions (e.g., CASI-5 domain severities)
- `Y`: non-specific observable behaviors (e.g., ABC items),

and learns a small set of orthogonal Aberrant Behavior Profiles (ABPs) -- linear combinations of items in `Y` --- that maximize sparsity of the differential severity --- profile associations. Concretely, after within-group rank standardization, DSV forms the group-differential association matrix between `X` and `Y`, projects it onto a low-dimensional orthonormal basis of `Y`, and then applies a varimax rotation to concentrate differential signal onto minimally overlapping subsets of severity dimensions. This construction yields ABPs that remain interpretable at the item level while supporting principled hypothesis testing of group-differential associations.

# Installation

> library(devtools)

> install_github("ericstrobl/DSV")

> library(DSV)

# Inputs

- `X` : n x pX matrix of severity scores / predictors (rows = subjects)

- `Y` : n x pY matrix of behaviors / items (rows = subjects)

- `labels` : length-n vector of two group labels (e.g., "ASD" and "TD")
Used for within-group rank-standardization and for computing differential associations.

- `verbal` : length-n vector (kept for interface consistency; may be used in extensions)

- `nc` : number of retained components (dimensions) to extract

- `ee` (optional) : eigen(cov(Y)) object to reuse the same PCA basis across refits/permutations

- `MR_ref` (optional) : reference loading matrix used for orthogonal Procrustes alignment
(stabilizes sign/permutation ambiguity across refits, especially in permutation testing)



