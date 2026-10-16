library(jsonlite)
source("analysis/glicko2.R")

# per-wave runs use 75 sessions per period.
# thankfully found a glicko2 algo someone else had done in R 
# and just needed to make a few adjustments :)

comparisons_path <- "analysis/output/comparisons.csv"
manifest_path <- "stimuli/faces.json"
output_path <- "analysis/output/face_glicko2_scores.csv"

sessions_per_period <- 100000   # one period, since the stimuli are static

# canonical glicko2 widens a competitor's deviation when they sit out a period,
# because people get rustier. faces don't, so that won't be necessary here
inflate_inactive <- FALSE

TAU <- 0.5
EPSILON <- 1e-6
INITIAL_RATING <- 1500
INITIAL_RD <- 350
INITIAL_VOLATILITY <- 0.06
ALGORITHM_VERSION <- "glicko2-static-stimulus-v1" # needed to make adapting this easier 
# still replicates reasonably well

as_logical_safe <- function(x) {
  if (is.logical(x)) return(x)
  text <- tolower(trimws(as.character(x)))
  out <- rep(NA, length(text))
  out[text %in% c("true", "t", "1", "yes", "y")] <- TRUE
  out[text %in% c("false", "f", "0", "no", "n", "")] <- FALSE
  out
}


# read in data
comparisons <- read.csv(comparisons_path, stringsAsFactors = FALSE,
                        check.names = FALSE,
                        na.strings = c("", "NA", "null", "NULL"))
manifest <- fromJSON(manifest_path, simplifyDataFrame = TRUE)

comparisons$judgment_target <- tolower(trimws(comparisons$judgment_target))
comparisons$left_score <- suppressWarnings(as.numeric(comparisons$left_score))
comparisons$right_score <- suppressWarnings(as.numeric(comparisons$right_score))

for (column in c("test_mode", "general_glicko2_update_eligible",
                 "within_female_glicko2_update_eligible",
                 "within_male_glicko2_update_eligible")) {
  comparisons[[column]] <- as_logical_safe(comparisons[[column]])
}


# comparions that can move ratings
usable <- comparisons[
  comparisons$task == "main_comparison" &
  comparisons$judgment_target %in% c("masculine", "feminine") &
  !ifelse(is.na(comparisons$test_mode), FALSE, comparisons$test_mode) &
  (is.na(comparisons$completion_status) | comparisons$completion_status == "complete"), ]

# the two sides must carry a real result that adds to one
usable$score_values_valid <- ifelse(is.na(
  usable$left_score %in% c(0, 0.5, 1) &
  usable$right_score %in% c(0, 0.5, 1) &
  abs((usable$left_score + usable$right_score) - 1) < 1e-10), FALSE,
  usable$left_score %in% c(0, 0.5, 1) &
  usable$right_score %in% c(0, 0.5, 1) &
  abs((usable$left_score + usable$right_score) - 1) < 1e-10)


# order sessions by when they ran, number them within each question asked, then
# cut them into blocks of sessions_per_period
when <- ifelse(!is.na(usable$session_start) & nzchar(usable$session_start),
               as.character(usable$session_start),
               as.character(usable$record_created_at))

sessions <- unique(data.frame(session_id = usable$session_id,
                              judgment_target = usable$judgment_target,
                              started = when, stringsAsFactors = FALSE))
sessions$source_order <- seq_len(nrow(sessions))
sessions <- sessions[order(sessions$judgment_target,
                           is.na(sessions$started),
                           sessions$started,
                           sessions$source_order), ]
sessions$index_within_target <- ave(sessions$source_order, sessions$judgment_target,
                                    FUN = seq_along)
sessions$rating_period <- ceiling(sessions$index_within_target / sessions_per_period)

usable <- merge(usable, sessions[c("session_id", "judgment_target", "rating_period")],
                by = c("session_id", "judgment_target"), all.x = TRUE, sort = FALSE)


# which faces belong to which systems
manifest$face_id <- as.character(manifest$face_id)
manifest$predicted_gender <- as.integer(manifest$predicted_gender)
manifest$active <- as_logical_safe(manifest$active)
manifest$practice_only <- as_logical_safe(manifest$practice_only)

in_play <- manifest[manifest$active %in% TRUE & !(manifest$practice_only %in% TRUE), ]

scopes <- list(
  general = list(
    faces = in_play$face_id,
    pair_types = c("female_within", "male_within", "cross_gender"),
    eligible = "general_glicko2_update_eligible"
  ),
  within_female_predicted = list(
    faces = in_play$face_id[in_play$predicted_gender == 0],
    pair_types = "female_within",
    eligible = "within_female_glicko2_update_eligible"
  ),
  within_male_predicted = list(
    faces = in_play$face_id[in_play$predicted_gender == 1],
    pair_types = "male_within",
    eligible = "within_male_glicko2_update_eligible"
  )
)


# presentations and response counts come from every trial the face appeared in;
# wins, losses and ties only from the trials that may move a rating
face_summary <- function(all_trials, rating_trials, faces) {
  rows <- lapply(faces, function(face_id) {
    shown_left <- all_trials[all_trials$left_face_id == face_id, ]
    shown_right <- all_trials[all_trials$right_face_id == face_id, ]
    rated_left <- rating_trials[rating_trials$left_face_id == face_id, ]
    rated_right <- rating_trials[rating_trials$right_face_id == face_id, ]

    scores <- c(rated_left$left_score, rated_right$right_score)
    opponents <- c(rated_left$right_face_id, rated_right$left_face_id)
    responses <- c(shown_left$response_type, shown_right$response_type)

    data.frame(
      face_id = face_id,
      n_presentations = nrow(shown_left) + nrow(shown_right),
      n_rating_eligible = length(scores),
      n_unique_opponents = length(unique(opponents[!is.na(opponents)])),
      n_wins = sum(scores == 1, na.rm = TRUE),
      n_losses = sum(scores == 0, na.rm = TRUE),
      n_ties = sum(scores == 0.5, na.rm = TRUE),
      n_unsure = sum(responses == "unsure", na.rm = TRUE),
      n_missing_or_invalid = sum(!(responses %in% c("left", "right", "equal", "unsure")),
                                 na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}


# run 
all_scores <- list()
all_history <- list()

for (trait in c("masculine", "feminine")) {
  for (scope_name in names(scopes)) {

    scope <- scopes[[scope_name]]

    all_trials <- usable[usable$judgment_target == trait &
                         usable$pair_type %in% scope$pair_types &
                         usable$left_face_id %in% scope$faces &
                         usable$right_face_id %in% scope$faces, ]

    rating_trials <- all_trials[all_trials[[scope$eligible]] %in% TRUE &
                                all_trials$score_values_valid, ]

    state <- data.frame(face_id = as.character(scope$faces),
                        rating = INITIAL_RATING,
                        rating_deviation = INITIAL_RD,
                        volatility = INITIAL_VOLATILITY,
                        stringsAsFactors = FALSE)

    periods <- sort(unique(usable$rating_period[usable$judgment_target == trait]))
    periods <- periods[is.finite(periods)]
    if (length(periods) == 0) periods <- 1

    for (period in periods) {
      games <- rating_trials[rating_trials$rating_period == period,
                             c("left_face_id", "right_face_id",
                               "left_score", "right_score")]
      state <- glicko2_update_period(state, games, tau = TAU, epsilon = EPSILON,
                                     inflate_inactive = inflate_inactive)

      snapshot <- state
      snapshot$judgment_target <- trait
      snapshot$rating_scope <- scope_name
      snapshot$rating_period <- period
      snapshot$n_games_in_period <- nrow(games)
      all_history[[length(all_history) + 1]] <- snapshot
    }

    result <- merge(state, face_summary(all_trials, rating_trials, scope$faces),
                    by = "face_id", all.x = TRUE, sort = FALSE)
    result$judgment_target <- trait
    result$rating_scope <- scope_name

    # a face nobody rated has no rank, rather than a rank based on its starting value
    rated <- !is.na(result$n_rating_eligible) & result$n_rating_eligible > 0
    result$rating_rank <- NA_integer_
    result$rating_rank[rated] <- rank(-result$rating[rated], ties.method = "min")

    result$n_rating_periods <- length(periods)
    result$sessions_per_period_per_trait <- sessions_per_period
    result$tau <- TAU
    result$initial_rating <- INITIAL_RATING
    result$initial_rating_deviation <- INITIAL_RD
    result$initial_volatility <- INITIAL_VOLATILITY
    result$inactive_rd_inflation <- inflate_inactive
    result$rating_algorithm_version <- ALGORITHM_VERSION

    all_scores[[length(all_scores) + 1]] <- result
  }
}

scores <- do.call(rbind, all_scores)
history <- do.call(rbind, all_history)

# carry the face metadata through, so the scores table stands on its own
scores <- merge(scores, in_play[setdiff(names(in_play), c("active", "practice_only"))],
                by = "face_id", all.x = TRUE, sort = FALSE)
scores <- scores[order(scores$judgment_target, scores$rating_scope, scores$rating_rank), ]

write.csv(scores, output_path, row.names = FALSE, na = "")
write.csv(history, sub("\\.csv$", "_history.csv", output_path),
          row.names = FALSE, na = "")
