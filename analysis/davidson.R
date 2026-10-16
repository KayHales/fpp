# thankfully this was used by other people, so i just had to adapt a few things :) 
# validated it with football outcomes and pokemon TCG
SCALE <- 173.7178 # rating points per logit
START <- 1500 


# everything the likelihood needs for one set of comparisons.
# parameters are packed as: strengths, side effect, log nu, rater tie offsets
davidson_setup <- function(data, faces, raters) {

  left_of <- match(data$left_face_id, faces)
  right_of <- match(data$right_face_id, faces)
  rater_of <- if (length(raters) > 0) match(data$session_id, raters) else NA
  y <- data$left_score
  is_tie <- y == 0.5

  n_face <- length(faces)
  n_rater <- length(raters)

  at_side <- n_face + 1L
  at_tie <- n_face + 2L
  at_offsets <- if (n_rater > 0) (n_face + 3L):(n_face + 2L + n_rater) else integer(0)
  n_par <- n_face + 2L + n_rater

  probabilities <- function(par) {
    a <- par[left_of] + par[at_side]
    b <- par[right_of]
    log_nu <- par[at_tie] + if (n_rater > 0) par[at_offsets][rater_of] else 0

    exp_a <- exp(a)
    exp_b <- exp(b)
    exp_tie <- exp(log_nu + (a + b) / 2)
    total <- exp_a + exp_b + exp_tie

    list(left = exp_a / total, right = exp_b / total, tie = exp_tie / total)
  }

  log_likelihood <- function(par) {
    p <- probabilities(par)
    sum(ifelse(is_tie, log(pmax(p$tie, 1e-300)),
        ifelse(y == 1, log(pmax(p$left, 1e-300)), log(pmax(p$right, 1e-300)))))
  }

  # negative log likelihood plus the ridge penalty
  objective <- function(par, penalty_strength, penalty_offset) {
    penalty <- 0.5 * penalty_strength * sum(par[seq_len(n_face)]^2)
    if (n_rater > 0) {
      penalty <- penalty + 0.5 * penalty_offset * sum(par[at_offsets]^2)
    }
    -(log_likelihood(par) - penalty)
  }

  gradient <- function(par, penalty_strength, penalty_offset) {
    p <- probabilities(par)
    score_left <- y - (p$left + 0.5 * p$tie)
    score_right <- (1 - y) - (p$right + 0.5 * p$tie)
    score_tie <- as.numeric(is_tie) - p$tie

    g <- numeric(n_par)
    g[seq_len(n_face)] <-
      as.vector(tapply(score_left, factor(left_of, levels = seq_len(n_face)),
                       sum, default = 0)) +
      as.vector(tapply(score_right, factor(right_of, levels = seq_len(n_face)),
                       sum, default = 0))
    g[at_side] <- sum(score_left)
    g[at_tie] <- sum(score_tie)
    if (n_rater > 0) {
      g[at_offsets] <- as.vector(tapply(score_tie,
                                        factor(rater_of, levels = seq_len(n_rater)),
                                        sum, default = 0))
      g[at_offsets] <- g[at_offsets] - penalty_offset * par[at_offsets]
    }
    g[seq_len(n_face)] <- g[seq_len(n_face)] - penalty_strength * par[seq_len(n_face)]
    -g
  }

  # observed information, inverted later to get standard errors
  hessian <- function(par, penalty_strength, penalty_offset) {
    p <- probabilities(par)
    pl <- p$left; pr <- p$right; pt <- p$tie

    d_aa <- pl * (1 - pl - 0.5 * pt) + 0.5 * pt * (0.5 - pl - 0.5 * pt)
    d_bb <- pr * (1 - pr - 0.5 * pt) + 0.5 * pt * (0.5 - pr - 0.5 * pt)
    d_ab <- -pl * (pr + 0.5 * pt) + 0.5 * pt * (0.5 - pr - 0.5 * pt)
    d_tt <- pt * (1 - pt)
    d_at <- -pl * pt + 0.5 * pt * (1 - pt)
    d_bt <- -pr * pt + 0.5 * pt * (1 - pt)

    H <- matrix(0, n_par, n_par)
    for (row in seq_len(nrow(data))) {
      i <- left_of[row]
      j <- right_of[row]

      H[i, i] <- H[i, i] + d_aa[row]
      H[j, j] <- H[j, j] + d_bb[row]
      H[i, j] <- H[i, j] + d_ab[row]
      H[j, i] <- H[j, i] + d_ab[row]

      H[i, at_side] <- H[i, at_side] + d_aa[row]
      H[at_side, i] <- H[at_side, i] + d_aa[row]
      H[j, at_side] <- H[j, at_side] + d_ab[row]
      H[at_side, j] <- H[at_side, j] + d_ab[row]
      H[at_side, at_side] <- H[at_side, at_side] + d_aa[row]

      H[i, at_tie] <- H[i, at_tie] + d_at[row]
      H[at_tie, i] <- H[at_tie, i] + d_at[row]
      H[j, at_tie] <- H[j, at_tie] + d_bt[row]
      H[at_tie, j] <- H[at_tie, j] + d_bt[row]
      H[at_tie, at_tie] <- H[at_tie, at_tie] + d_tt[row]
    }

    diag(H)[seq_len(n_face)] <- diag(H)[seq_len(n_face)] + penalty_strength
    if (n_rater > 0) diag(H)[at_offsets] <- diag(H)[at_offsets] + penalty_offset
    H
  }

  list(n_par = n_par, n_face = n_face, n_rater = n_rater,
       at_side = at_side, at_tie = at_tie, at_offsets = at_offsets,
       log_likelihood = log_likelihood, objective = objective,
       gradient = gradient, hessian = hessian)
}


davidson_fit <- function(data, faces, raters, penalty_strength,
                         penalty_offset = 1) {
  model <- davidson_setup(data, faces, raters)
  found <- optim(numeric(model$n_par),
                 fn = function(p) model$objective(p, penalty_strength, penalty_offset),
                 gr = function(p) model$gradient(p, penalty_strength, penalty_offset),
                 method = "BFGS",
                 control = list(maxit = 600, reltol = 1e-11))
  list(par = found$par, model = model, converged = found$convergence == 0)
}


# pick the ridge penalty by held-out likelihood, splitting raters rather than
# rows. a held-out rater is scored with the global tie propensity
davidson_choose_penalty <- function(data, faces,
                                    grid = c(0.1, 0.25, 0.5, 1, 2, 4, 8),
                                    n_folds = 5, seed = 20260924) {
  raters <- sort(unique(data$session_id))
  set.seed(seed)
  fold_of <- sample(rep(seq_len(n_folds), length.out = length(raters)))
  names(fold_of) <- raters

  held_out <- numeric(length(grid))

  for (g in seq_along(grid)) {
    total <- 0
    for (fold in seq_len(n_folds)) {
      train_raters <- raters[fold_of[raters] != fold]
      train <- data[data$session_id %in% train_raters, ]
      test <- data[!data$session_id %in% train_raters, ]
      if (nrow(train) == 0 || nrow(test) == 0) next

      fitted <- davidson_fit(train, faces, train_raters, grid[g], 1)

      scoring <- davidson_setup(test, faces, character(0))
      par <- numeric(scoring$n_par)
      par[seq_along(faces)] <- fitted$par[seq_along(faces)]
      par[scoring$at_side] <- fitted$par[fitted$model$at_side]
      par[scoring$at_tie] <- fitted$par[fitted$model$at_tie]
      total <- total + scoring$log_likelihood(par)
    }
    held_out[g] <- total
  }

  grid[which.max(held_out)]
}
