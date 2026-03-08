
library(tidyverse)
library(lmerTest)
library(mice)
library(metafor)
library(marginaleffects)
library(broom.mixed)
current_path <- getwd()
message("当前执行路径：", current_path)
processed_data_dir <- file.path(current_path, "data", "processed")
external_factors_path <- file.path(current_path, "data", "country_level_factors.xlsx")
list.files(
  path = processed_data_dir, 
  pattern = "\\.RData$",
  full.names = TRUE
) %>% 
  walk(~load(.x, envir = .GlobalEnv)) 
regional_metrics <- read_excel(external_factors_path)
print(summary(regional_metrics))



region_codes <- c(40, 56, 76, 156, 191, 203, 208, 233, 250, 276, 300, 376, 380, 442, 484, 528, 705, 724, 752, 756, 840)
region_names <- c("Austria", "Belgium", "Brazil", "China", "Croatia", "Czech Republic", "Denmark", "Estonia", "France", "Germany", "Greece", "Israel", "Italy", "Luxembourg", "Mexico", "Netherlands", "Slovenia", "Spain", "Sweden", "Switzerland", "USA")


process_mi_list <- function(idx, mi_lists) {
  combined <- bind_rows(
    mi_lists$hrs[[idx]], 
    mi_lists$charls[[idx]], 
    mi_lists$share[[idx]], 
    mi_lists$mhas[[idx]]
  ) %>%
    mutate(
      id = as.character(id),

      country_label = factor(country, levels = region_codes, labels = region_names)
    )
  return(combined)
}
chess_mi_merged <- map(1:10, ~process_mi_list(.x, list_of_mi_sources))

run_chess_lmm <- function(data_list, country_name, outcome_var, adjusted = FALSE) {
  base_formula <- paste(outcome_var, "~ chess_habit * cwave + (1 | id)")
  adj_covariates <- "age + male + lowedu + lowwealth + lbrf + nonmarried + hhres + contact + smoken + weekdrink + phyinact + dis + adl_any + iadl"
  formula_str <- if(adjusted) paste(base_formula, "+", adj_covariates) else base_formula
  target_data_list <- map(data_list, ~filter(.x, country_label == country_name))
  first_df <- target_data_list[[1]]
  if (all(!is.na(first_df))) {
    model <- lmer(as.formula(formula_str), data = first_df, control = lmerControl(optimizer = "bobyqa"))
    res <- avg_comparisons(model, variables = "chess_habit", re.form = NA) %>% tidy()
  } else {
    models <- map(target_data_list, ~lmer(as.formula(formula_str), data = .x, control = lmerControl(optimizer = "bobyqa")))
    res <- map(models, ~avg_comparisons(.x, variables = "chess_habit", re.form = NA)) %>%
      mice::pool() %>%
      tidy()
  }
  return(res %>% mutate(country = country_name))
}

cog_results_adj <- map_dfr(region_names, ~run_chess_lmm(chess_mi_merged, .x, "cog_score", adjusted = TRUE))
meta_cog <- rma(yi = estimate, sei = std.error, data = cog_results_adj, slab = country)

meta_reg_data <- cog_results_adj %>%
  left_join(country_factor_data, by = "country")

res_mod_skills <- rma(yi = estimate, sei = std.error, 
                      mods = ~ poly(digital_skills, 2), 
                      data = meta_reg_data)

ggplot(cog_results_adj, aes(x = estimate, y = reorder(country, estimate))) +
  geom_point(size = 3) +
  geom_errorbarh(aes(xmin = estimate - 1.96*std.error, xmax = estimate + 1.96*std.error), height = .2) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  theme_minimal() +
  labs(title = "Effect of Chess Habit on Cognitive Score", x = "Estimate (95% CI)", y = "")


chess_forest_theme <- forestploter::forest_theme(
  core = list(bg_params = list(fill = c("white"))), 
  base_size = 10,
  ci_pch = 18,              
  ci_col = "#2A52BE",      
  ci_alpha = 0.9,
  ci_lty = 1,
  ci_lwd = 1.2,
  ci_Theight = unit(0, "cm"),
  xaxis_lwd = 1,
  xaxis_cex = 0.9,
  summary_col = "#C04000", 
  refline_col = "darkgray", 
  refline_lwd = 1.5, 
  title_just = "center",
  title_cex = 1.1,
  title_fontface = "bold",
  arrow_type = "closed",
  arrow_label_just = "end",
  arrow_length = 0.04,
  arrow_fill = "gray20"
)


prepare_forest_data <- function(results_df, meta_model) {
  dat <- results_df %>%
    mutate(weights = weights(meta_model)) %>%
    bind_rows(
      tibble(country = "Summary (Random Effects)", estimate = as.numeric(meta_model$beta),
             std.error = as.numeric(meta_model$se), weights = NA),
      tibble(country = NA), tibble(country = NA), tibble(country = NA)
    ) %>%
    mutate(
      lower = estimate - 1.96 * std.error,
      upper = estimate + 1.96 * std.error,
      weights_text = ifelse(is.na(weights), "", sprintf("%0.2f", weights)),
      beta_ci = ifelse(is.na(estimate), "", sprintf("%0.2f [%0.2f, %0.2f]", estimate, lower, upper)),
      blank_space = paste(rep(" ", 20), collapse = " "),
      plot_slot = paste(rep(" ", 30), collapse = " ")
    )
  
  fig_df <- dat %>% select(country, blank_space, plot_slot, beta_ci, weights_text)
  fig_df[is.na(fig_df)] <- " "
  colnames(fig_df) <- c("Region", "", "Effect Size", "AME (95% CI)", "Weight (%)")
  
  return(list(full = dat, fig = fig_df))
}

cog_data <- prepare_forest_data(cog_results_adj, meta_cog_model)

p_cog <- forestploter::forest(
  cog_data$fig,
  est = cog_data$full$estimate,
  lower = cog_data$full$lower, 
  upper = cog_data$full$upper,
  sizes = c(sqrt(cog_data$full$weights[1:(nrow(cog_data$full)-4)]/5), rep(1.2, 4)),
  is_summary = c(rep(FALSE, nrow(cog_data$full)-4), TRUE, rep(FALSE, 3)),
  ci_column = 3,
  ref_line = 0,
  arrow_lab = c("Lower cognitive score", "Higher cognitive score"), 
  xlim = c(-0.4, 0.4),
  theme = chess_forest_theme
)


add_stat_annotations <- function(plot_obj, meta_obj, row_idx) {
  txt_het <- bquote(paste("Heterogeneity: ", italic(tau)^2, "=", .(sprintf("%.2f", meta_obj$tau2)), 
                          "; ", italic(I)^2, "=", .(sprintf("%.1f", meta_obj$I2)), "%"))
  
  txt_overall <- bquote(paste("Overall effect: ", italic(z), "=", .(sprintf("%.2f", meta_obj$zval)), 
                              ", ", italic(P), .(ifelse(meta_obj$pval < 0.001, " < 0.001", sprintf(" = %.3f", meta_obj$pval)))))
  
  plot_obj <- add_text(plot_obj, text = txt_het, row = row_idx - 2, col = 1, just = "left", gp = gpar(fontsize = 8), parse = TRUE)
  plot_obj <- add_text(plot_obj, text = txt_overall, row = row_idx, col = 1, just = "left", gp = gpar(fontsize = 8), parse = TRUE)
  return(plot_obj)
}

p_cog <- add_stat_annotations(p_cog, meta_cog_model, nrow(cog_data$full))

combined_forest <- wrap_elements(p_cog) + wrap_elements(p_satlife) + wrap_elements(p_shlt) + 
  plot_layout(ncol = 3) + 
  plot_annotation(tag_levels = "A", title = "Association between Chess Habits and Health Outcomes")



