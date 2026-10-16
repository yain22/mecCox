method_map <- make_method_label_map(res = res, Lin_Wei_ps_method = Lin_Wei_ps_method, include_unweighted_cox = include_unweighted_cox, include_mec_no_cf = include_mec_no_cf)

method_levels <- method_map$Method

method_short_levels <- method_map$Method_Short

sample_grid_plot <- summ %>% distinct(n1, n0) %>% mutate(N_total = n1 + n0) %>% arrange(N_total)

plot_dat <- tidyr::expand_grid(sample_grid_plot, Method = factor(method_levels, levels = method_levels)) %>% left_join(summ %>% mutate(N_total = n1 + n0, Method = factor(Method, levels = method_levels)), by = c("n1", "n0", "N_total", "Method")) %>% left_join(method_map, by = "Method") %>% mutate(Method_Short = factor(Method_Short, levels = method_short_levels)) %>% arrange(Method_Short, N_total)

