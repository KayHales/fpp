source("analysis/davidson.R")

comparisons <- read.csv("analysis/output/comparisons.csv", stringsAsFactors = FALSE)

comparisons <- comparisons[comparisons$completion_status == "complete", ]
comparisons$left_score <- suppressWarnings(as.numeric(comparisons$left_score))
comparisons <- comparisons[comparisons$left_score %in% c(0, 0.5, 1), ]

scopes <- list(
  general = c("female_within", "male_within", "cross_gender"),
  within_female_predicted = "female_within",
  within_male_predicted = "male_within"
)

scores <- data.frame()
settings <- data.frame()

for (this_target in c("masculine", "feminine")) {
  for (scope_name in names(scopes)) {

    subset <- comparisons[comparisons$judgment_target == this_target &
                          comparisons$pair_type %in% scopes[[scope_name]], ]
    if (nrow(subset) == 0) next

    faces <- sort(unique(c(subset$left_face_id, subset$right_face_id)))
    raters <- sort(unique(subset$session_id))

    penalty <- davidson_choose_penalty(subset, faces)
    fitted <- davidson_fit(subset, faces, raters, penalty)

    strength <- fitted$par[seq_along(faces)]

    # standard errors from the inverse of the observed information
    H <- fitted$model$hessian(fitted$par, penalty, 1)
    covariance <- try(solve(H), silent = TRUE)
    standard_error <- rep(NA_real_, length(faces))
    if (!inherits(covariance, "try-error")) {
      variance <- diag(covariance)[seq_along(faces)]
      standard_error[variance > 0] <- sqrt(variance[variance > 0])
    }

    appearances <- as.vector(table(factor(c(subset$left_face_id, subset$right_face_id),
                                         levels = faces)))

    scores <- rbind(scores, data.frame(
      face_id = faces,
      judgment_target = this_target,
      rating_scope = scope_name,
      rating = START + SCALE * strength,
      rating_se = SCALE * standard_error,
      strength_logit = strength,
      n_comparisons = appearances,
      stringsAsFactors = FALSE
    ))

    settings <- rbind(settings, data.frame(
      judgment_target = this_target,
      rating_scope = scope_name,
      n_comparisons = nrow(subset),
      n_faces = length(faces),
      lambda_strength = penalty,
      side_effect_logit = fitted$par[fitted$model$at_side],
      side_effect_rating_pts = SCALE * fitted$par[fitted$model$at_side],
      tie_nu_global = exp(fitted$par[fitted$model$at_tie]),
      logLik = fitted$model$log_likelihood(fitted$par),
      converged = fitted$converged,
      stringsAsFactors = FALSE
    ))
  }
}

# rank within each system, best first
scores$rating_rank <- NA_integer_
for (key in unique(paste(scores$judgment_target, scores$rating_scope))) {
  rows <- paste(scores$judgment_target, scores$rating_scope) == key
  scores$rating_rank[rows] <- rank(-scores$rating[rows], ties.method = "min")
}

write.csv(scores, "analysis/output/face_davidson_scores.csv", row.names = FALSE, na = "")
write.csv(settings, "analysis/output/face_davidson_scores_model.csv", row.names = FALSE, na = "")
