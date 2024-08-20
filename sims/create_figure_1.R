# create Figure 1 (example VIM trajectories) and compute true VIM values of all summaries

library("dplyr")
library("broom")
library("purrr")
library("tidyr")
library("ggplot2")
library("cowplot")
library("lvimp")
theme_set(theme_cowplot())

# create the true VIMs ---------------------------------------------------------
vims_1 <- tibble::tibble(s = "1", k = 1:5,
                         vim = c(1:5 * 0.1))
vims_2 <- tibble::tibble(s = "2", k = 1:5,
                         vim = rep(0.8, 5))
vims_3 <- tibble::tibble(s = "3", k = 1:5,
                         vim = c(2:0 * 0.1, rep(0, 2)))
all_vims <- tibble::as_tibble(rbind.data.frame(vims_1, vims_2, vims_3))

# create the plot --------------------------------------------------------------
longitudinal_vim_plot <- all_vims %>%
  ggplot(aes(x = k, y = vim, shape = s)) +
  geom_point(size = 3) +
  geom_line(linetype = "dashed") +
  labs(x = "Time point", y = "True VIM", shape = "Feature of interest") +
  theme(legend.position = "bottom", legend.direction = "horizontal")
ggsave(here::here("..", "..", "plots", "lvim_example.png"),
       longitudinal_vim_plot, width = 4, height = 3, units = "in", dpi = 300)

# true VIM summaries -----------------------------------------------------------
# mean over the entire time series
means <- all_vims |> 
  group_by(s) |> 
  summarize(vim_summary = mean(vim)) |> 
  mutate(summary_type = "mean", .before = vim_summary)
Orange %>% split(.$Tree) %>% map(~lm(age ~ 1 + circumference, data = .x)) %>% map_df(tidy) %>% filter(term == 'circumference')
trends <- all_vims %>%
  split(.$s) |> 
  map(~lm(vim ~ k, data = .x)) |> 
  map_df(tidy) |> 
  select(term, estimate) |>
  mutate(s = c("1", "1", "2", "2", "3", "3")) |> 
  mutate(summary_type = ifelse(term == "(Intercept)", "trend_intercept", "trend_slope")) |> 
  select(s, summary_type, estimate) |> 
  rename(vim_summary = estimate)
  

piecewise_linear_estimate <- function(x) {
  if (!is.matrix(x)) {
    indices <- seq_len(length(x))
    x <- matrix(x, nrow = 1)
  }
  indices <- seq_len(ncol(x))
  return(x[, range(indices, na.rm = TRUE)[1]] / 2 +
           x[, range(indices, na.rm = TRUE)[2]] / 2 +
           rowSums(x[, 2:(range(indices, na.rm = TRUE)[2] - 1), drop = FALSE], na.rm = TRUE))
}
gm_mean <- function(x, na.rm = TRUE, zero.propagate = FALSE){
  if (any(x < 0, na.rm = TRUE)){
    return(NaN)
  }
  if (zero.propagate){
    if(any(x == 0, na.rm = TRUE)){
      return(0)
    }
    exp(mean(log(x), na.rm = na.rm))
  } else {
    exp(sum(log(x[x > 0]), na.rm = na.rm) / length(x))
  }
}
get_gmrc <- function(x, deriv_func) {
  abs_deriv <- abs(deriv_func(x))
  if (!any(abs_deriv > 0)) {
    return(0)
  } else {
    return(exp(mean(log(abs_deriv), na.rm = TRUE)))
  }
}
# get_gmrc <- function(x, deriv_func, ...) {
#   return(gm_mean(deriv_func(x), ...))
# }
autc_piecewise_linear <- tibble::tibble(
  "s" = c("1", "2", "3"),
  "summary_type" = "autc_piecewise_linear",
  "vim_summary" = c(
    piecewise_linear_estimate(all_vims$vim[all_vims$s == "1"]),
    piecewise_linear_estimate(all_vims$vim[all_vims$s == "2"]),
    piecewise_linear_estimate(all_vims$vim[all_vims$s == "3"])
  )
)
spline_1 <- splinefun(x = all_vims$k[all_vims$s == "1"], 
                      y = all_vims$vim[all_vims$s == "1"])
spline_2 <- splinefun(x = all_vims$k[all_vims$s == "2"], 
                      y = all_vims$vim[all_vims$s == "2"])
spline_3 <- splinefun(x = all_vims$k[all_vims$s == "3"], 
                      y = all_vims$vim[all_vims$s == "3"])
autc_spline <- tibble::tibble(
  "s" = c("1", "2", "3"),
  "summary_type" = "autc_spline",
  "vim_summary" = c(
    integrate(spline_1,
              lower = range(all_vims$k[all_vims$s == "1"])[1],
              upper = range(all_vims$k[all_vims$s == "1"])[2])$value,
    integrate(spline_2,
              lower = range(all_vims$k[all_vims$s == "2"])[1],
              upper = range(all_vims$k[all_vims$s == "2"])[2])$value,
    integrate(spline_3,
              lower = range(all_vims$k[all_vims$s == "3"])[1],
              upper = range(all_vims$k[all_vims$s == "3"])[2])$value
  )
)
piecewise_linear_deriv <- function(x, vims) {
  return(ifelse(x < length(vims), vims[x + 1] - vims[x], 0))
}
grid_len <- 1000
grid <- seq(range(all_vims$k[all_vims$s == "1"])[1], range(all_vims$k[all_vims$s == "1"])[2],
            length.out = grid_len)
gmrc_piecewise_linear <- tibble::tibble(
  "s" = c("1", "2", "3"),
  "summary_type" = "gmrc_piecewise_linear",
  "vim_summary" = c(
    get_gmrc(grid[-grid_len], deriv_func = function(x) piecewise_linear_deriv(x, vims = all_vims$vim[all_vims$s == "1"])),
    get_gmrc(grid[-grid_len], deriv_func = function(x) piecewise_linear_deriv(x, vims = all_vims$vim[all_vims$s == "2"])),
    get_gmrc(grid[-grid_len], deriv_func = function(x) piecewise_linear_deriv(x, vims = all_vims$vim[all_vims$s == "3"]))
  )
)

gmrc_spline <- tibble::tibble(
  "s" = c("1", "2", "3"),
  "summary_type" = "gmrc_spline",
  "vim_summary" = c(
    get_gmrc(grid[-grid_len], deriv_func = function(x) spline_1(x, deriv = 1)),
    get_gmrc(grid[-grid_len], deriv_func = function(x) spline_2(x, deriv = 1)),
    get_gmrc(grid[-grid_len], deriv_func = function(x) spline_3(x, deriv = 1))
  )
)

rbind(means, trends, autc_piecewise_linear, autc_spline, gmrc_piecewise_linear, gmrc_spline) |> 
  print(n = Inf)
