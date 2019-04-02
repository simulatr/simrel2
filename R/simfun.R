# ## Essential Function ----
# bdiag <- function(a_mat, b_mat) {
#     da   <- dim(a_mat)
#     db   <- dim(b_mat)
#     up   <- cbind(a_mat, matrix(0, da[1], db[2]))
#     down <- cbind(matrix(0, db[1], da[2]), b_mat)
#     return(rbind(up, down))
# }
# compute_eigen <- function(i, decay) {
#     exp(-1 * decay * (i - 1))
# }
# rotation_block <- function(n) {
#     Qmat <- matrix(rnorm(n ^ 2), n)
#     Qmat <- scale(Qmat, scale = FALSE)
#     qr.Q(qr(Qmat))
# }
# pred_pos <- function(relpos, p, q) {
#     if (depth(relpos) == 1) {
#         relpos <- list(relpos)
#     }
#     n_relpred     <- q - sapply(relpos, length)
#     n_irrelpred   <- p - sum(q)
#     n_fctr        <- rep(seq.int(n_relpred), n_relpred)
#     names(relpos) <- seq.int(relpos)
#     extra_pos     <- setdiff(1:p, unlist(relpos))
#     new_pos       <- split(sample(extra_pos, sum(n_relpred)), n_fctr)
#     irrel_pos     <- setdiff(extra_pos, Reduce(union, new_pos))
#     relpred       <- lapply(names(relpos), function(x) c(relpos[[x]], new_pos[[x]]))
#     out           <- append(relpred, list(irrel_pos))
#     Map(sort, Filter(length, out))
# }
# rot_mat <- function(pred_pos_list) {
#     idx <- order(unlist(pred_pos_list))
#     n   <- sapply(pred_pos_list, length)
#     out <- lapply(n, rotation_block)
#     Reduce(bdiag, out)[idx, idx]
# }
# get_cov <- function(pos, Rsq, eta, p, lambda){
#     out      <- vector("numeric", p)
#     alph     <- runif(length(pos), -1, 1)
#     out[pos] <- sign(alph) * sqrt(Rsq * abs(alph) / sum(abs(alph)) * lambda[pos] * eta)
#     return(out)
# }
# cov_mat <- function(relpos, rsq, p, m, eta, gamma) {
#     if (depth(relpos) == 1) {
#         relpos <- list(relpos)
#     }
#     lambda    <- sapply(1:p, compute_eigen, decay = gamma)
#     kappa     <- sapply(1:m, compute_eigen, decay = eta)
#     rel_cov   <- sapply(seq_along(relpos), function(idx) {
#         get_cov(relpos[[idx]], rsq[idx], kappa[idx], p, lambda)
#     })
#     irrel_cov <- matrix(0, nrow = p, ncol = m - length(relpos))
#     cbind(rel_cov, irrel_cov)
# }
# get_sigma <- function(sigma_ww, cov_zw, sigma_zz) {
#     rbind(
#         cbind(sigma_ww, t(cov_zw)),
#         cbind(cov_zw, sigma_zz)
#     )
# }
#
# ## Extra Properties Functions ----
# beta_z <- function(lambda, cov_xy) {
#     diag(1/lambda) %*% cov_xy
# }
# beta <- function(rot_x, beta_z, rot_y) {
#     Reduce(`%*%`, list(rot_x, beta_z, rot_y))
# }
# beta0 <- function(beta, mu_x = NULL, mu_y = NULL) {
#     result <- rep(0, ncol(beta))
#     if (!is.null(mu_y)) result <- result + mu_y
#     if (!is.null(mu_x)) result <- result - t(beta) %*% mu_x
# }
# rsq_w <- function(cov_zw, lambda, kappa) {
#     var_w      <- diag(1/sqrt(kappa))
#     sigma_zinv <- diag(1/lambda)
#     Reduce(`%*%`, list(var_w, t(cov_zw), sigma_zinv, cov_zw, var_w))
# }
# rsq_y <- function(cov_zw, lambda, kappa, rot_y) {
#     sigma_yy   <- Reduce(`%*%`, list(t(rot_y), diag(kappa), rot_y))
#     sigma_zinv <- diag(1/lambda)
#     var_y      <- diag(1/sqrt(diag(sigma_yy)))
#     Reduce(`%*%`, list(var_y, rot_y, t(cov_zw), sigma_zinv, cov_zw, t(rot_y), var_y))
# }
# minerror <- function(rot_y, cov_zw, lambda, kappa) {
#     sigma_ww   <- diag(kappa)
#     sigma_zinv <- diag(1/lambda)
#     expl_var   <- Reduce(`%*%`, list(t(cov_zw), sigma_zinv, cov_zw))
#     Reduce(`%*%`, list(t(rot_y), (sigma_ww - expl_var), rot_y))
# }
# get_data <- function(n, p, m, sigma, rot_x, rot_y = NULL, mu_x = NULL, mu_y = NULL){
#     rotate_y  <- all(exists("rot_y"), !is.null(rot_y))
#     sigma_rot <- chol(sigma)
#     train_cal <- matrix(rnorm(n * (p + m), 0, 1), nrow = n) %*% sigma_rot
#     Z         <- train_cal[, (m + 1):(m + p), drop = F]
#     X         <- Z %*% t(rot_x)
#     W         <- train_cal[, 1:m, drop = F]
#     Y         <- if (rotate_y) W %*% t(rot_y) else W
#
#     if (!is.null(mu_x)) X <- sweep(X, 2, mu_x, "+")
#     if (!is.null(mu_y)) Y <- sweep(Y, 2, mu_y, "+")
#
#     list(y = unname(Y), x = X)
# }
#
# ## Wrapper Function ----
# simrel <- function(p, q, m, gamma, eta, relpos, ypos, rsq,
#                    mu_x = NULL, mu_y = NULL, seed = 123) {
#     set.seed(seed)
#     original_seed <- seed
#     rm(seed, envir = environment())
#
#     ## Parse the lists
#     relpos  <- chr2list(relpos)
#     ypos    <- chr2list(ypos)
#     q       <- chr2list(q)
#     rsq     <- chr2list(rsq)
#
#     ## Fix some null values for uni-variate simulation
#     if (is.null(eta)) eta <- 1
#     if (is.null(ypos)) ypos <- matrix(1, 1, 1)
#
#     ## Eigenvalues of response and predictors
#     eigen_x <- compute_eigen(1:p, gamma)
#     eigen_y <- compute_eigen(1:m, eta)
#
#     ## Sample and fix relevant predictors (randomness involved)
#     relpred <- pred_pos(relpos, p, q)
#
#     ## Construct covariance matrix (randomness involved)
#     cov_zw <- cov_mat(relpos, rsq, p, m, eta, gamma)
#
#     ## Construct variance-covariance matrices
#     var_ww   <- diag(eigen_y)
#     var_zz   <- diag(eigen_x)
#     var_zinv <- diag(1/eigen_x)
#     sigma_zw <- get_sigma(var_ww, cov_zw, var_zz)
#
#     ## Rotation matrix from standard normal (randomness involved)
#     rot_x <- as.matrix(rot_mat(relpred))
#     rot_y <- as.matrix(rot_mat(ypos))
#
#     ## Variance and covariance of response and predictors
#     var_yy <- Reduce(`%*%`, list(rot_y, var_ww, t(rot_y)))
#     var_xx <- Reduce(`%*%`, list(rot_x, var_zz, t(rot_x)))
#     cov_xy <- Reduce(`%*%`, list(rot_x, cov_zw, t(rot_y)))
#
#     ## Extra methods for various properties
#     sigma  <- function() get_sigma(var_yy, cov_xy, var_xx)
#     beta_z <- function() var_zinv %*% cov_zw
#     beta   <- function() Reduce(`%*%`, list(rot_x, beta_z(), t(rot_y)))
#     beta0  <- function(mu_x = NULL, mu_y = NULL) {
#         beta <- beta()
#         result <- as.matrix(rep(0, ncol(beta)))
#         if (!is.null(mu_y)) result <- result + mu_y
#         if (!is.null(mu_x)) result <- result - t(beta) %*% mu_x
#         return(result)
#     }
#     minerror <- function() {
#         expl_var <- Reduce(`%*%`, list(t(cov_zw), var_zinv, cov_zw))
#         Reduce(`%*%`, list(t(rot_y), (var_ww - expl_var), rot_y))
#     }
#     data <- function(n, mu_x = NULL, mu_y = NULL, seed = original_seed) {
#         set.seed(seed)
#         get_data(n, p, m, sigma_zw, rot_x, rot_y, mu_x = mu_x, mu_y = mu_y)
#     }
#     rsq_w <- function() {
#         var_w <- diag(1/sqrt(eigen_y))
#         Reduce(`%*%`, list(var_w, t(cov_zw), var_zinv, cov_zw, var_w))
#     }
#     rsq_y <- function() {
#         sigma_yy <- Reduce(`%*%`, list(t(rot_y), diag(eigen_y), rot_y))
#         var_y_elem <- -1/sqrt(diag(sigma_yy))
#         var_y    <- diag(var_y_elem, length(var_y_elem), length(var_y_elem))
#         Reduce(`%*%`, list(var_y, rot_y, t(cov_zw), var_zinv, cov_zw, t(rot_y), var_y))
#     }
#     out <- list(
#         relpred      = relpred,
#         sigma        = sigma,
#         rotation_x   = rot_x,
#         rotation_y   = rot_y,
#         beta         = beta,
#         beta0        = beta0,
#         minerror     = minerror,
#         r_squared_w  = rsq_w,
#         r_squared    = rsq_y,
#         data         = data
#     )
#     class(out) <- "simrel"
#     return(out)
# }
#
# ## Methods ------------------
# print.simrel <- function(obj, ...) {
#     print(ls.str(obj))
# }
#
# ## Implementation ---------------------
# sobj_m <- simrel(
#     p = 15,
#     q = "3, 4, 5",
#     m = 5,
#     gamma = 0.7,
#     eta = 0.4,
#     relpos = "1:3; 4:6; 7:10",
#     ypos = "1, 5; 2, 4; 3",
#     rsq = "0.4, 0.6, 0.8",
#     seed = 123
# )
# sobj_1 <- simrel(
#     p = 10,
#     q = "5",
#     m = 1,
#     gamma = 0.7,
#     eta = NULL,
#     relpos = "1:5",
#     ypos = NULL,
#     rsq = "0.8",
#     seed = 123
# )
#
# cat("\nMulti-Response Linear Model\n")
# print(sobj_m)
# cat("\nUni-Response Linear Model\n")
# print(sobj_1)
