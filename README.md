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

# Inputs to DSV

- `X` : n x pX matrix of severity scores / predictors (rows = subjects)

- `Y` : n x pY matrix of behaviors / items (rows = subjects)

- `labels` : length-n vector of two group labels (e.g., "ASD" and "TD")
Used for within-group rank-standardization and for computing differential associations.

- `verbal` : length-n vector

- `nc` : number of retained components (dimensions) to extract

- `ee` (optional) : eigen(cov(Y)) object to reuse the same PCA basis across refits/permutations

- `MR_ref` (optional) : reference loading matrix used for orthogonal Procrustes alignment
(stabilizes sign/permutation ambiguity across refits, especially in permutation testing)

# Run the Algorithm

> n  <- 200 # synthetic data

> pX <- 12

> pY <- 58

> X <- matrix(rnorm(n * pX), n, pX)

> Y <- matrix(rnorm(n * pY), n, pY)

> colnames(X) <- paste0("CASI_", seq_len(pX))

> colnames(Y) <- paste0("ABC_",  seq_len(pY))

> labels <- rep(c("ASD", "TD"), length.out = n) # two-group labels (e.g., ASD vs. TD)

> pa <- parallel_analysis_DSV(Y, labels, B = 200, alpha = 0.95) # parallel analysis

> nc <- pa$k

> nc

> mod <- DSV(X, Y, labels, nc = nc) # run DSV algorithm

> dim(mod$MR)       # pX x nc

> dim(mod$RtW)      # nc x pY

> dim(mod$optimal_components)  # n x nc

> res <- permutation_testing_DSV(X, Y, labels, nperms = 1000, nc = nc) # permutation testing

> res$factors$stat

> res$factors$pval

> res$factors$qB

# Outputs

DSV returns a list containing (most commonly used items first):

- `MR` : pX x nc matrix
Differential severity–component associations
(Spearman correlation in group 1 minus Spearman correlation in group 2, after rank-standardization)

- `RtW` : nc x pY matrix
Rotated component–item weight map (ABP definitions).
You can interpret each row as an item-weighted ABP.

- `optimal_components` : n x nc matrix
Subject-level rotated component scores.

- `effects` : pX x pY matrix
Implied differential effects from severity dimensions to items through the component map (MR %*% RtW).

- `R` : nc x nc orthogonal varimax rotation

- `Q` : nc x nc orthogonal Procrustes alignment (identity if MR_ref is NULL)

- `eigen` : eigen-decomposition used for the PCA basis of Y

- `X`, `Y` : rank-standardized versions used internally

- `total_variance_Y` : trace of cov(Y) (used for variance accounting)

- `loadings_unrot`, `loadings_rot` : pY x nc loading-style matrices before/after rotation

- `Vaccounted_unrot`, `Vaccounted_rot` : variance accounting summaries


permutation_testing_DSV returns a list containing:

- `omnibus` : global statistic and p-value across all MR entries

- `factors` : component-wise stats, permutation p-values, maxT FWER p-values, and FDR q-values

- `MR_entrywise` : entrywise |MR| stats, permutation p-values, FDR q-values, and FDR curve
