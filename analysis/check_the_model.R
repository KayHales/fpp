source("analysis/davidson.R")

# fit on four fifths of the raters, predict the other fifth, and compare the
# predicted probabilities to what those raters really chose. 

comparisons <- read.csv("analysis/output/comparisons.csv", stringsAsFactors = FALSE)

comparisons <- comparisons[comparisons$completion_status == "complete", ]
comparisons$left_score <- suppressWarnings(as.numeric(comparisons$left_score))
comparisons <- comparisons[comparisons$left_score %in% c(0, 0.5, 1), ]
settings <- read.csv("analysis/output/face_davidson_scores_model.csv", stringsAsFactors = FALSE)

n_folds <- 5
seed <- 20261014

scopes <- list(
  general = c("female_within", "male_within", "cross_gender"),
  within_female_predicted = "female_within",
  within_male_predicted = "male_within"
)

# 95% CIs
wilson <- function(hits, n) {
  if (n == 0) return(c(NA_real_, NA_real_))
  p <- hits / n
  z <- 1.959964
  denom <- 1 + z^2 / n
  centre <- (p + z^2 / (2 * n)) / denom
  halfwidth <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / denom
  c(centre - halfwidth, centre + halfwidth)
}

bins <- data.frame()
summaries <- data.frame()

for (this_target in c("masculine", "feminine")) {
  for (scope_name in names(scopes)) {

    subset <- comparisons[comparisons$judgment_target == this_target &
                          comparisons$pair_type %in% scopes[[scope_name]], ]
    if (nrow(subset) == 0) next

    faces <- sort(unique(c(subset$left_face_id, subset$right_face_id)))
    raters <- sort(unique(subset$session_id))
    penalty <- settings$lambda_strength[settings$judgment_target == this_target &
                                settings$rating_scope == scope_name]

    set.seed(seed)
    fold_of <- sample(rep(seq_len(n_folds), length.out = length(raters)))
    names(fold_of) <- raters

    p_left <- rep(NA_real_, nrow(subset))
    p_right <- rep(NA_real_, nrow(subset))
    p_tie <- rep(NA_real_, nrow(subset))
    gap <- rep(NA_real_, nrow(subset))

    for (fold in seq_len(n_folds)) {
      train_raters <- raters[fold_of[raters] != fold]
      train <- subset[subset$session_id %in% train_raters, ]
      test_rows <- which(!subset$session_id %in% train_raters)
      if (length(test_rows) == 0 || nrow(train) == 0) next

      fitted <- davidson_fit(train, faces, train_raters, penalty)
      strength <- fitted$par[seq_along(faces)]
      side_effect <- fitted$par[fitted$model$at_side]
      log_nu <- fitted$par[fitted$model$at_tie]
      train_offsets <- fitted$par[fitted$model$at_offsets]

      a <- strength[match(subset$left_face_id[test_rows], faces)] + side_effect
      b <- strength[match(subset$right_face_id[test_rows], faces)]
      exp_a <- exp(a)
      exp_b <- exp(b)
      half <- exp((a + b) / 2)
      sum_left <- 0; sum_right <- 0; sum_tie <- 0
      for (offset in train_offsets) {
        exp_tie <- exp(log_nu + offset) * half
        total <- exp_a + exp_b + exp_tie
        sum_left <- sum_left + exp_a / total
        sum_right <- sum_right + exp_b / total
        sum_tie <- sum_tie + exp_tie / total
      }
      n_offsets <- length(train_offsets)

      p_left[test_rows] <- sum_left / n_offsets
      p_right[test_rows] <- sum_right / n_offsets
      p_tie[test_rows] <- sum_tie / n_offsets
      gap[test_rows] <- strength[match(subset$left_face_id[test_rows], faces)] -
                        strength[match(subset$right_face_id[test_rows], faces)]
    }

    have <- !is.na(p_left)
    y <- subset$left_score[have]
    p_left <- p_left[have]; p_right <- p_right[have]
    p_tie <- p_tie[have]; gap <- gap[have]

    add_bin <- function(check, bin, n, predicted, observed, hits, gap_mid = NA) {
      interval <- wilson(hits, n)
      data.frame(judgment_target = this_target, rating_scope = scope_name, check = check,
                 bin = bin, n = n, predicted = predicted, observed = observed,
                 lo = interval[1], hi = interval[2], gap_mid = gap_mid,
                 stringsAsFactors = FALSE)
    }

    # 1. among comparisons somebody actually decided, does the predicted chance
    # the left face wins match how often it won?
    decided <- y != 0.5
    p_decided <- p_left[decided] / (p_left[decided] + p_right[decided])
    won <- as.numeric(y[decided] == 1)
    which_bin <- cut(p_decided, quantile(p_decided, seq(0, 1, 0.1)),
                     include.lowest = TRUE, labels = FALSE)
    for (k in sort(unique(which_bin))) {
      rows <- which_bin == k
      bins <- rbind(bins, add_bin("decisive", k, sum(rows),
                                  mean(p_decided[rows]), mean(won[rows]), sum(won[rows])))
    }

    # 2. does the predicted tie probability match how often people tied?
    which_bin <- cut(p_tie, quantile(p_tie, seq(0, 1, 0.1)),
                     include.lowest = TRUE, labels = FALSE)
    for (k in sort(unique(which_bin))) {
      rows <- which_bin == k
      bins <- rbind(bins, add_bin("tie", k, sum(rows), mean(p_tie[rows]),
                                  mean(y[rows] == 0.5), sum(y[rows] == 0.5)))
    }

    # 3. is the model right about when ties happen?
    distance <- abs(gap)
    which_bin <- cut(distance, quantile(distance, seq(0, 1, length.out = 9)),
                     include.lowest = TRUE, labels = FALSE)
    for (k in sort(unique(which_bin))) {
      rows <- which_bin == k
      bins <- rbind(bins, add_bin("tie_by_gap", k, sum(rows), mean(p_tie[rows]),
                                  mean(y[rows] == 0.5), sum(y[rows] == 0.5),
                                  mean(distance[rows]) * SCALE))
    }

    # regress what happened on the logit of what was predicted. 
    safe <- pmin(pmax(p_decided, 1e-6), 1 - 1e-6)
    logit_predicted <- log(safe / (1 - safe))
    recalibration <- glm(won ~ logit_predicted, family = binomial())

    # and how much better than knowing only the overall rates of each outcome
    base_rates <- c(mean(y == 1), mean(y == 0.5), mean(y == 0))
    model_fit <- mean(ifelse(y == 0.5, log(p_tie),
                      ifelse(y == 1, log(p_left), log(p_right))))
    base_fit <- mean(ifelse(y == 0.5, log(base_rates[2]),
                     ifelse(y == 1, log(base_rates[1]), log(base_rates[3]))))

    summaries <- rbind(summaries, data.frame(
      judgment_target = this_target,
      rating_scope = scope_name,
      n_comparisons = length(y),
      n_faces = length(faces),
      penalty = penalty,
      calibration_slope = unname(coef(recalibration)[2]),
      calibration_intercept = unname(coef(recalibration)[1]),
      tie_rate_observed = mean(y == 0.5),
      tie_rate_predicted = mean(p_tie),
      logloss_model = -model_fit,
      logloss_baserate = -base_fit,
      stringsAsFactors = FALSE
    ))
  }
}

write.csv(bins, "analysis/output/calibration_bins.csv", row.names = FALSE)
write.csv(summaries, "analysis/output/calibration.csv", row.names = FALSE)
summaries
