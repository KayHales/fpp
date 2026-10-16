source("analysis/davidson.R")

# how much would the scores change if we had asked a different set of people?

comparisons <- read.csv("analysis/output/comparisons.csv", stringsAsFactors = FALSE)

comparisons <- comparisons[comparisons$completion_status == "complete", ]
comparisons$left_score <- suppressWarnings(as.numeric(comparisons$left_score))
comparisons <- comparisons[comparisons$left_score %in% c(0, 0.5, 1), ]
settings <- read.csv("analysis/output/face_davidson_scores_model.csv", stringsAsFactors = FALSE)

n_splits <- 6
project_to <- 230   # this can change, it's just for projecting to various sample sizes
seed <- 20260925

scopes <- list(
  general = c("female_within", "male_within", "cross_gender"),
  within_female_predicted = "female_within",
  within_male_predicted = "male_within"
)

set.seed(seed)
splits <- data.frame()

for (split in seq_len(n_splits)) {
  for (this_target in sort(unique(comparisons$judgment_target))) {

    for_target <- comparisons[comparisons$judgment_target == this_target, ]
    # sorted so the split depends on the seed and not on the row order of the file
    sessions <- sort(unique(for_target$session_id))
    session_wave <- for_target$wave[match(sessions, for_target$session_id)]

    # take half of each wave, so both halves cover every cell
    first_half <- c()
    for (w in sort(unique(session_wave))) {
      in_wave <- sessions[session_wave == w]
      first_half <- c(first_half, sample(in_wave, floor(length(in_wave) / 2)))
    }

    for (scope_name in names(scopes)) {

      subset <- for_target[for_target$pair_type %in% scopes[[scope_name]], ]
      faces <- sort(unique(c(subset$left_face_id, subset$right_face_id)))

      penalty <- settings$lambda_strength[settings$judgment_target == this_target &
                                  settings$rating_scope == scope_name]
      if (length(penalty) != 1) penalty <- 0.5

      a <- subset[subset$session_id %in% first_half, ]
      b <- subset[!subset$session_id %in% first_half, ]

      strength_a <- davidson_fit(a, faces, sort(unique(a$session_id)),
                                 penalty)$par[seq_along(faces)]
      strength_b <- davidson_fit(b, faces, sort(unique(b$session_id)),
                                 penalty)$par[seq_along(faces)]

      half_r <- cor(strength_a, strength_b, method = "spearman")

      splits <- rbind(splits, data.frame(
        split = split,
        judgment_target = this_target,
        rating_scope = scope_name,
        model = "Davidson-BT",
        n_raters = length(sessions),
        reliability = 2 * half_r / (1 + half_r),
        stringsAsFactors = FALSE
      ))
    }
  }
}


# average over the splits, and keep the range as a sense of how much it moved
summary_stats <- aggregate(reliability ~ judgment_target + rating_scope + model + n_raters, splits,
                           function(x) c(mean = mean(x), lo = min(x), hi = max(x)))
results <- data.frame(summary_stats[, 1:4],
                      reliability = summary_stats$reliability[, "mean"],
                      lo = summary_stats$reliability[, "lo"],
                      hi = summary_stats$reliability[, "hi"])

# and project to a different rater count, same formula, k times as many raters
k <- project_to / results$n_raters
results$project_to_raters <- project_to
results$reliability_projected <- k * results$reliability / (1 + (k - 1) * results$reliability)
results$n_splits <- n_splits

results <- results[order(-results$reliability), ]

write.csv(results, "analysis/output/reliability.csv", row.names = FALSE)
results
