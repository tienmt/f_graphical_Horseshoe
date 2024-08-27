# tempered graphical Horsshoe
simulations for tempered/fractional graphical Horsshoe. See paper: "Mai, T. T. (2024). Concentration of a sparse Bayesian model with Horseshoe prior in estimating high-dimensional precision matrix. arXiv:2406.14269."

* need to load file 'Untitled.R' and 'Omega0_p100_rnd.Rdata'.
* Then run the file 'sim_dense.R' for simulations with dense precision matrices.
* run the file 'simu_sparse.R' for simulations with sparse precision matrices.

The main Gibbs sampler is due to "Li, Y., Craig, B. A., & Bhadra, A. (2019). The graphical horseshoe estimator for inverse covariance matrices. Journal of Computational
and Graphical Statistics, 28(3), 747–757."

The sparse precision matrix in the file 'Omega0_p100_rnd.Rdata' is from the paper "Zhang, R., Yao, Y., & Ghosh, M. (2022). Contraction of a quasi-bayesian model with shrinkage priors in precision matrix estimation. Journal of Statistical Planning and Inference, 221, 154–171."

Anindya Bhadra is the author of "The graphical horseshoe estimator for inverse covariance matrices. JCGS" as well as "Sagar, K., Banerjee, S., Datta, J., & Bhadra, A. (2024). Precision matrix estimation under the horseshoe-like prior–penalty dual. EJS" and the package 'GHS'. It was nice for me to meet and talk with him at ISBA 2024 in Venice.

