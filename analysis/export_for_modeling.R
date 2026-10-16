# make the csv for the ML model
pooled_dir <- "analysis/output"
OUT <- "modeling_handoff"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

d <- read.csv(file.path(pooled_dir, "face_davidson_scores.csv"), stringsAsFactors = FALSE)
g <- read.csv(file.path(pooled_dir, "face_glicko2_scores.csv"), stringsAsFactors = FALSE)
rel <- read.csv(file.path(pooled_dir, "reliability.csv"), stringsAsFactors = FALSE)
cmp <- read.csv(file.path(pooled_dir, "comparisons.csv"), stringsAsFactors = FALSE)
cmp <- cmp[cmp$task == "main_comparison" & cmp$completion_status == "complete", ]

O <- pooled_dir  


# data
meta <- unique(g[, c("face_id", "filename", "image_url", "predicted_gender",
                     "predicted_gender_label", "reliability_score_gender",
                     "skin_tone", "skin_tone_group",
                     "predicted_age", "target", "intersection")])
names(meta)[names(meta) == "target"] <- "is_ai_generated"

wide <- meta

for (tr in c("masculine", "feminine")) {
  for (sc in c("general", "within_female_predicted", "within_male_predicted")) {
    x <- d[d$judgment_target == tr & d$rating_scope == sc,
           c("face_id", "rating", "rating_se", "n_comparisons")]
    stub <- paste0(substr(tr, 1, 4), "_", sub("_predicted", "", sc))
    names(x) <- c("face_id", paste0(stub, c("", "_se", "_n")))
    wide <- merge(wide, x, by = "face_id", all.x = TRUE)
  }
}

# masculinity and femininity are one bipolar axis, so the reversed average of
# the two standardised scales is a single more reliable target per face
for (sc in c("general", "within_female_predicted", "within_male_predicted")) {
  m <- d[d$judgment_target == "masculine" & d$rating_scope == sc, c("face_id", "rating")]
  f <- d[d$judgment_target == "feminine" & d$rating_scope == sc, c("face_id", "rating")]
  z <- merge(m, f, by = "face_id", suffixes = c("_m", "_f"))
  z$comp <- (as.vector(scale(z$rating_m)) - as.vector(scale(z$rating_f))) / 2
  z <- z[, c("face_id", "comp")]
  names(z)[2] <- paste0("composite_", sub("_predicted", "", sc))
  wide <- merge(wide, z, by = "face_id", all.x = TRUE)
}

# the 20 anchor faces are never rated, so they have no score in any system and
# must not become training rows
wide$has_scores <- !is.na(wide$masc_general)

write.csv(wide, file.path(OUT, "face_targets.csv"), row.names = FALSE, na = "")

long <- merge(d[, c("face_id", "judgment_target", "rating_scope", "rating",
                    "rating_se", "rating_rank", "n_comparisons")],
              meta, by = "face_id")
write.csv(long, file.path(OUT, "face_targets_long.csv"), row.names = FALSE, na = "")
