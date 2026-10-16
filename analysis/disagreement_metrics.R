# how much do raters disagree about each face?

comparisons_path <- "analysis/output/comparisons.csv"
output_path <- "analysis/output/face_disagreement_metrics.csv"

comparisons <- read.csv(comparisons_path, stringsAsFactors = FALSE)

as_logical_safe <- function(x) {
  if (is.logical(x)) return(x)
  text <- tolower(trimws(as.character(x)))
  out <- rep(NA, length(text))
  out[text %in% c("true", "t", "1", "yes", "y")] <- TRUE
  out[text %in% c("false", "f", "0", "no", "n")] <- FALSE
  out
}

entropy <- function(counts) {
  total <- sum(counts, na.rm = TRUE)
  if (total <= 0) return(NA_real_)
  p <- counts[counts > 0] / total
  -sum(p * log2(p))
}

comparisons$test_mode <- as_logical_safe(comparisons$test_mode)

main <- comparisons[
  comparisons$task == "main_comparison" &
  comparisons$judgment_target %in% c("masculine", "feminine") &
  !ifelse(is.na(comparisons$test_mode), FALSE, comparisons$test_mode) &
  (is.na(comparisons$completion_status) | comparisons$completion_status == "complete"), ]

scopes <- list(
  general = c("female_within", "male_within", "cross_gender"),
  within_female_predicted = "female_within",
  within_male_predicted = "male_within"
)


# one face's record within one system.
# a face appears on both sides of the screen, so "won" means chosen when it was
# on the left and chosen when it was on the right.
summarise_face <- function(data, face_id, trait, scope) {

  as_left <- data[data$left_face_id == face_id, ]
  as_right <- data[data$right_face_id == face_id, ]

  outcomes <- c(
    ifelse(as_left$response_type == "left", "win",
    ifelse(as_left$response_type == "right", "loss",
    ifelse(as_left$response_type == "equal", "tie",
    ifelse(as_left$response_type == "unsure", "unsure", "invalid")))),
    ifelse(as_right$response_type == "right", "win",
    ifelse(as_right$response_type == "left", "loss",
    ifelse(as_right$response_type == "equal", "tie",
    ifelse(as_right$response_type == "unsure", "unsure", "invalid"))))
  )

  opponents <- c(as_left$right_face_id, as_right$left_face_id)

  n_win <- sum(outcomes == "win", na.rm = TRUE)
  n_loss <- sum(outcomes == "loss", na.rm = TRUE)
  n_tie <- sum(outcomes == "tie", na.rm = TRUE)
  n_unsure <- sum(outcomes == "unsure", na.rm = TRUE)
  n_total <- length(outcomes)
  n_decisive <- n_win + n_loss

  data.frame(
    face_id = face_id,
    judgment_target = trait,
    rating_scope = scope,
    n_presentations = n_total,
    n_unique_opponents = length(unique(opponents[!is.na(opponents)])),
    n_wins = n_win,
    n_losses = n_loss,
    n_ties = n_tie,
    n_unsure = n_unsure,
    n_invalid = sum(outcomes == "invalid", na.rm = TRUE),
    decisive_win_proportion = if (n_decisive > 0) n_win / n_decisive else NA_real_,
    tie_rate = if (n_total > 0) n_tie / n_total else NA_real_,
    unsure_rate = if (n_total > 0) n_unsure / n_total else NA_real_,
    binary_decisive_entropy_bits = entropy(c(n_win, n_loss)),
    substantive_response_entropy_bits = entropy(c(n_win, n_loss, n_tie)),
    all_response_entropy_bits = entropy(c(n_win, n_loss, n_tie, n_unsure)),
    maximum_substantive_entropy_bits = if (n_win + n_loss + n_tie > 0) log2(3) else NA_real_,
    stringsAsFactors = FALSE
  )
}


rows <- list()

for (trait in c("masculine", "feminine")) {
  for (scope_name in names(scopes)) {

    subset <- main[main$judgment_target == trait &
                   main$pair_type %in% scopes[[scope_name]], ]

    faces <- sort(unique(c(subset$left_face_id, subset$right_face_id)))
    faces <- faces[!is.na(faces) & nzchar(faces)]

    for (face_id in faces) {
      rows[[length(rows) + 1]] <- summarise_face(subset, face_id, trait, scope_name)
    }
  }
}

write.csv(do.call(rbind, rows), output_path, row.names = FALSE, na = "")
