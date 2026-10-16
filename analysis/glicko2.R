# adapted for our study :) see the score_glicko2.R and the davidson.R script
GLICKO2_SCALE <- 173.7178


# how much an opponent's uncertainty discounts the information from beating them
g_of <- function(phi) 1 / sqrt(1 + 3 * phi^2 / pi^2)

# chance of beating an opponent at mu_j, given their uncertainty
expected_score <- function(mu, mu_j, phi_j) {
  1 / (1 + exp(-g_of(phi_j) * (mu - mu_j)))
}

solve_volatility <- function(phi, sigma, delta, v, tau, epsilon = 1e-6) {

  a <- log(sigma^2)

  f <- function(x) {
    ex <- exp(x)
    ex * (delta^2 - phi^2 - v - ex) / (2 * (phi^2 + v + ex)^2) - (x - a) / tau^2
  }

  # bracket the root between A and B
  A <- a
  if (delta^2 > phi^2 + v) {
    B <- log(delta^2 - phi^2 - v)
  } else {
    # walk downwards in steps of tau until f turns positive
    k <- 1
    B <- a - k * tau
    while (f(B) < 0) {
      k <- k + 1
      B <- a - k * tau
    }
  }

  f_A <- f(A)
  f_B <- f(B)

  while (abs(B - A) > epsilon) {
    C <- A + (A - B) * f_A / (f_B - f_A)
    f_C <- f(C)

    if (f_C * f_B <= 0) {
      A <- B
      f_A <- f_B
    } else {
      f_A <- f_A / 2
    }

    B <- C
    f_B <- f_C
  }

  exp(A / 2)
}


# one face's update for one rating period.
glicko2_update_one <- function(rating, rating_deviation, volatility, opponents,
                               tau = 0.5, epsilon = 1e-6,
                               inflate_inactive = TRUE) {

  mu <- (rating - 1500) / GLICKO2_SCALE
  phi <- rating_deviation / GLICKO2_SCALE

  if (is.null(opponents) || nrow(opponents) == 0) {
    phi_new <- if (isTRUE(inflate_inactive)) sqrt(phi^2 + volatility^2) else phi
    return(c(rating = rating,
             rating_deviation = GLICKO2_SCALE * phi_new,
             volatility = volatility))
  }

  mu_j <- (opponents$rating - 1500) / GLICKO2_SCALE
  phi_j <- opponents$rating_deviation / GLICKO2_SCALE
  g_j <- g_of(phi_j)

  expected <- 1 / (1 + exp(-g_j * (mu - mu_j)))
  expected <- pmin(pmax(expected, 1e-15), 1 - 1e-15)

  v <- 1 / sum(g_j^2 * expected * (1 - expected))

  delta <- v * sum(g_j * (opponents$score - expected))

  volatility_new <- solve_volatility(phi, volatility, delta, v, tau, epsilon)

  phi_star <- sqrt(phi^2 + volatility_new^2)
  phi_new <- 1 / sqrt(1 / phi_star^2 + 1 / v)
  mu_new <- mu + phi_new^2 * sum(g_j * (opponents$score - expected))

  c(rating = GLICKO2_SCALE * mu_new + 1500,
    rating_deviation = GLICKO2_SCALE * phi_new,
    volatility = volatility_new)
}

glicko2_update_period <- function(state, games, tau = 0.5, epsilon = 1e-6,
                                  inflate_inactive = TRUE) {

  before <- state
  rownames(before) <- before$face_id
  after <- before

  for (face_id in before$face_id) {

    as_left <- games[games$left_face_id == face_id, ]
    as_right <- games[games$right_face_id == face_id, ]

    # who this face met, and how it did, from both sides of the pair
    met <- rbind(
      data.frame(opponent = as_left$right_face_id, score = as_left$left_score,
                 stringsAsFactors = FALSE),
      data.frame(opponent = as_right$left_face_id, score = as_right$right_score,
                 stringsAsFactors = FALSE)
    )

    opponents <- data.frame(
      rating = before[met$opponent, "rating"],
      rating_deviation = before[met$opponent, "rating_deviation"],
      score = met$score,
      stringsAsFactors = FALSE
    )

    updated <- glicko2_update_one(
      rating = before[face_id, "rating"],
      rating_deviation = before[face_id, "rating_deviation"],
      volatility = before[face_id, "volatility"],
      opponents = opponents,
      tau = tau, epsilon = epsilon, inflate_inactive = inflate_inactive
    )

    after[face_id, "rating"] <- updated[["rating"]]
    after[face_id, "rating_deviation"] <- updated[["rating_deviation"]]
    after[face_id, "volatility"] <- updated[["volatility"]]
  }

  rownames(after) <- NULL
  after
}
