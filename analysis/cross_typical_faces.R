input_path <- "analysis/output/face_glicko2_scores.csv"
output_path <- "analysis/output/cross_typical_faces.csv"
top_n <- 25

scores <- read.csv(input_path, stringsAsFactors = FALSE, check.names = FALSE,
                   na.strings = c("", "NA", "null", "NULL"))

scores$face_id <- as.character(scores$face_id)
scores$judgment_target <- tolower(trimws(scores$judgment_target))
scores$rating_scope <- tolower(trimws(scores$rating_scope))
for (column in c("predicted_gender", "rating", "rating_deviation", "volatility",
                 "rating_rank", "n_rating_eligible", "n_unique_opponents")) {
  scores[[column]] <- suppressWarnings(as.numeric(scores[[column]]))
}

# face metadata rides along on the Glicko output, so carry through whatever
# of it is present
metadata <- intersect(
  c("filename", "image_url", "predicted_gender_label", "predicted_age",
    "reliability_score_gender", "reliability_score_age", "skin_tone",
    "skin_tone_group", "skin_tone_group_label", "target", "target_label",
    "intersection", "intersection_label", "manifest_version",
    "rating_algorithm_version"),
  names(scores))

rating_columns <- c("rating", "rating_deviation", "volatility", "rating_rank",
                    "n_rating_eligible", "n_unique_opponents")


build_profile <- function(profile_name, trait, within_scope, gender_code, n_keep) {

  within <- scores[scores$judgment_target == trait &
                   scores$rating_scope == within_scope &
                   scores$predicted_gender == gender_code, ]

  general <- scores[scores$judgment_target == trait &
                    scores$rating_scope == "general" &
                    scores$predicted_gender == gender_code, ]

  # highest rated first; ties broken towards the better measured face, then the
  # more often rated one, then the id so the order is never arbitrary
  within <- within[order(-within$rating, within$rating_deviation,
                         -within$n_rating_eligible, within$face_id), ]
  within <- head(within, n_keep)

  general_lookup <- general[, c("face_id", rating_columns)]
  names(general_lookup)[-1] <- paste0("general_", names(general_lookup)[-1])

  out <- merge(within[, unique(c("face_id", metadata, rating_columns))],
               general_lookup, by = "face_id", all.x = TRUE, sort = FALSE)

  out <- out[order(match(out$face_id, within$face_id)), ]

  names(out)[match(rating_columns, names(out))] <- paste0("within_", rating_columns)

  out$profile <- profile_name
  out$trait <- trait
  out$within_scope <- within_scope
  out$profile_rank <- seq_len(nrow(out))
  out$approx_within_lower_95 <- out$within_rating - 1.96 * out$within_rating_deviation
  out$approx_within_upper_95 <- out$within_rating + 1.96 * out$within_rating_deviation

  leading <- c("profile", "profile_rank", "face_id", "trait", "within_scope",
               "within_rating", "within_rating_deviation", "within_volatility",
               "within_rating_rank", "approx_within_lower_95",
               "approx_within_upper_95", "within_n_rating_eligible",
               "within_n_unique_opponents", "general_rating",
               "general_rating_deviation", "general_volatility",
               "general_rating_rank", "general_n_rating_eligible",
               "general_n_unique_opponents")
  out[, c(leading, setdiff(names(out), leading))]
}


output <- rbind(
  build_profile("female_predicted_high_masculinity", "masculine",
                "within_female_predicted", 0, top_n),
  build_profile("male_predicted_high_femininity", "feminine",
                "within_male_predicted", 1, top_n)
)

write.csv(output, output_path, row.names = FALSE, na = "")
