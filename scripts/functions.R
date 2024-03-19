# Functions ------------------------------------------------------------------
myGplot.defaults = function(
    type = c("paper","poster","slides")[1],
    base_size = if (type == "paper") { 10 } else if (type == "slides") { 32 } else if (type == "poster") { 36 } else { 10 },
    margin=c("t" = 0.6, "r" = 0.5, "b" = 0.5, "l" = 0.3),
    set_theme = T
)
{
  require(ggplot2)
  t <- theme(axis.text.x = element_text(size=base_size, vjust=1),
             axis.text.y = element_text(size=base_size, hjust=1, vjust=.5),
             axis.title.x = element_text(size=base_size , vjust=0, hjust=0.5, face = "bold"),
             axis.title.y = element_text(size=base_size, hjust= 0.5, vjust=0.5, face = "bold"),
             strip.text = element_text(size=base_size, color = "white"),
             strip.background = element_rect(fill = "black", color = "black"),
             legend.title = element_text(size=base_size, face = "bold", hjust= 0),
             legend.text = element_text(size=base_size),
             plot.margin = unit(margin, "lines"),
             aspect.ratio = 1)
  
  if (set_theme) theme_set(theme_bw(base_size=base_size) + t) else return(t)
}

## General functions ------------------------------------------------------------------
se <- function(x) sqrt(var(x) / (length(x) - 1))

geometric.mean = psych::geometric.mean

get_vowel_stats <-
  . %>%
  select(-F0_gm) %>%
  summarise(
    across(ends_with("gm"),
           .fns = list("mean" = mean, "var" = var)),
    F1F2_cov = cov(F1_gm, F2_gm), 
    F1F3_cov = cov(F1_gm, F3_gm),
    F2F3_cov = cov(F2_gm, F3_gm),
    heightbackness_cov = cov(height_gm, backness_gm),
    synthtalkerF0_gm = synthtalker_F0)

## Normalization functions ------------------------------------------------------------
F0_to_SR = function(F0) {
  168 * (F0 / 168)^(1/3)
}

F1_to_height = function(F1, SR) {
  log(F1 / SR)
}

F2_to_backness = function(F1, F2) {
  log(F2 / F1)
}

height_to_F1 = function(height, SR) {
  exp(height) * SR
}

backness_to_F2 = function(backness, height, SR) {
  exp(backness) * height_to_F1(height, SR)
}

# Aliases
get_F1_Hz_from_height <- height_to_F1

get_F2_Hz_from_backness <- backness_to_F2

add_Hz_from_height_backness <- function(
    data,
    var.height = "Item.Height",
    var.backness = "Item.Backness",
    var.SR = "SR"
) {
  require(rlang)
  
  data %>%
    mutate(
      Item.Cue_Hz_F1 = height_to_F1(!! sym(var.height), !! sym(var.SR)),
      Item.Cue_Hz_F2 = backness_to_F2(!! sym(var.backness), !! sym(var.height), !! sym(var.SR)))
}

get_transformation <- function(
    data,
    cues = c("F0","F1", "F2", "F3")
) {
  data %>%
    mutate_at(
      cues,
      .funs =
        c("Mel" = function(x) phonR::normMel(x),
          "Bark" = function(x) phonR::normBark(x),
          "ERB" = function(x) phonR::normErb(x),
          "semitones" = function(x) 12 * log2(x/100)))
}

## Outlier detection and correction --------------------------------------------------
get_cumulative_probability = function(x1, x2, mean, sigma) {
  pmvnorm(lower = -Inf, upper = c(x1, x2), mean = mean, sigma = sigma)
}

get_cumulative_probability_allCues = function(x1, x2, x3, x4, x5, mean, sigma) {
  pmvnorm(lower = -Inf, upper = c(x1, x2, x3, x4, x5), mean = mean, sigma = sigma)
}

get_cumulative_probability_univariate = function(x, mean, sd) {
  pnorm(x, mean = mean, sd = sd)
}

is_outlier = function(x, cutoff = outlier_probability_cutoff) {
  !between(
    x, 
    cutoff / 2, 
    1 - cutoff / 2)
}

obtain_densities <- . %>%
  group_by(Talker, Vowel) %>%
  nest() %>%
  mutate(
    # density based on F1 and F2
    x_F1F2 = map(data, ~ cbind(.$F1, .$F2)),
    mean_F1F2 = map(x_F1F2, ~ colMeans(.x)),
    cov_F1F2 = map(x_F1F2, ~ cov(.x)),
    x_F1F2 = NULL) %>%
  unnest(data) %>%
  # normalize densities within each talker and vowel
  mutate(
    cumulative_probability_F1F2 = pmap(
      .l = list(F1, F2, mean_F1F2, cov_F1F2), 
      .f = get_cumulative_probability)) %>%
  mutate_at(
    vars(cumulative_probability_F1F2),
    unlist) %>%
  ungroup()

#Obtain univariate densities for outlier correction of all cues in raw Hz
obtain_densities_univariates <- . %>%
  group_by(Talker, Vowel) %>%
  nest() %>%
  mutate(
    # density based on F0
    x_F0 = map(data, ~ cbind(.$F0)),
    mean_F0 = map(x_F0, ~ mean(.x, na.rm = TRUE)),
    sd_F0 = map(x_F0, ~ sd(.x)),
    # density based on F1
    x_F1 = map(data, ~ cbind(.$F1)),
    mean_F1 = map(x_F1, ~ mean(.x)),
    sd_F1 = map(x_F1, ~ sd(.x)),
    # density based on F2
    x_F2 = map(data, ~ cbind(.$F2)),
    mean_F2 = map(x_F2, ~ mean(.x)),
    sd_F2 = map(x_F2, ~ sd(.x)),
    # density based on F3
    x_F3 = map(data, ~ cbind(.$F3)),
    mean_F3 = map(x_F3, ~ mean(.x)),
    sd_F3 = map(x_F3, ~ sd(.x)),
    # density based on Duration
    x_Duration = map(data, ~ cbind(.$Duration)),
    mean_Duration = map(x_Duration, ~ mean(.x)),
    sd_Duration = map(x_Duration, ~ sd(.x)),
    x_F0 = NULL,
    x_F1 = NULL,
    x_F2 = NULL,
    x_F3 = NULL,
    x_Duration = NULL) %>%
  unnest(data) %>%
  # normalize densities within each talker and vowel
  mutate(
    cumulative_probability_F0 = pmap(
      .l = list(F0, mean_F0, sd_F0), 
      .f = get_cumulative_probability_univariate),
    cumulative_probability_F1 = pmap(
      .l = list(F1, mean_F1, sd_F1), 
      .f = get_cumulative_probability_univariate),
    cumulative_probability_F2 = pmap(
      .l = list(F2, mean_F2, sd_F2),
      .f = get_cumulative_probability_univariate),
    cumulative_probability_F3 = pmap(
      .l = list(F3, mean_F3, sd_F3), 
      .f = get_cumulative_probability_univariate),
    cumulative_probability_F3 = pmap(
      .l = list(F3, mean_F3, sd_F3), 
      .f = get_cumulative_probability_univariate),
    cumulative_probability_Duration = pmap(
      .l = list(Duration, mean_Duration, sd_Duration), 
      .f = get_cumulative_probability_univariate)) %>%
  mutate_at(
    vars(starts_with("cumulative_probability")),
    unlist) %>%
  ungroup()

#Obtain multivariate densities for outlier correction of all cues in raw Hz
obtain_densities_allCues <- function(data) {
  data %>%
    group_by(Talker, Vowel) %>%
    nest() %>%
    mutate(
      # density based on all cues
      x_cues = map(data, ~ cbind(.$F0, .$F1, .$F2, .$F3, .$Duration)),
      mean_cues = map(x_cues, ~ colMeans(.x)),
      cov_cues = map(x_cues, ~ cov(.x))) %>%
    unnest(data) %>%
    # normalize densities within each talker and vowel
    mutate(
      cumulative_probability_allCues = pmap(
        .l = list(F0, F1, F2, F3, Duration, mean_cues, cov_cues), 
        .f = get_cumulative_probability_allCues)) %>%
    mutate_at(
      vars(cumulative_probability_allCues),
      unlist) %>%
    ungroup()
}

# Vowel plot function
plot_vowels <- function(data, x, y) {
  ggplot(data,
         aes(
           x = F2_gm, 
           y = F1_gm,
           color = Vowel)) +
    geom_point(alpha = .5) +
    scale_x_reverse() +
    scale_y_reverse() +
    vowelPlot_components
}
