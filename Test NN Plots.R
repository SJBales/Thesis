# MSE of Selection Criteria Across vars except gamma and beta
nn_analyzing_df %>%
  group_by(Method, Variable) %>%
  filter(Variable != 6 & Variable != 7, 
         norm_var == 'high') %>%
  summarize(MSE = mean(MSE)) %>%
  ggplot(aes(x = Variable, y = MSE, fill = Method)) +
  geom_bar(stat = 'identity', position = 'dodge')

nn_analyzing_df %>%
  filter(Variable == 6) %>%
  group_by(Method) %>%
  summarize(MSE = median(MSE)) %>%
  ggplot(aes(x = Method, y = MSE, fill = Method)) + geom_bar(stat = 'identity')

# MSE of criteria when the variability of normal predictors is high
nn_analyzing_df %>%
  group_by(Method, Variable) %>%
  filter(Variable != 6 & Variable != 7, 
         norm_var == 'high') %>%
  summarize(MSE = mean(MSE)) %>%
  ggplot(aes(x = Variable, y = MSE, fill = Method)) +
  geom_bar(stat = 'identity', position = 'dodge')

nn_analyzing_df %>%
  group_by(Method) %>%
  filter(Variable == 9) %>%
  summarize(MSE_mean = mean(MSE), 
            MSE_med = median(MSE),
            MSE_sd = sd(MSE),
            MSE_min = min(MSE),
            MSE_max = max(MSE))


library(ggridges)

# MSE Ridgeline for 
nn_analyzing_df %>%
  filter(Variable == 7) %>%
  ggplot(aes(x = MSE, y = Method, fill = Method)) +
  geom_density_ridges()

# Histograms of MSE for T dist
nn_analyzing_df %>%
  filter(Variable == 7, 
         norm_error == "low", 
         norm_var == "high") %>%
  ggplot(aes(x = MSE, fill = Method)) +
  geom_histogram() +
  facet_wrap(Method ~. )

nn_analyzing_df %>% 
  filter(Variable != 6 & Variable != 7) %>%
  group_by(Variable, norm_var, Method) %>%
  summarize(MSE = mean(MSE)) %>%
  ggplot(aes(x = Variable, y = MSE, fill = Method)) + 
  geom_bar(stat = "identity", position = 'dodge') + 
  facet_wrap(. ~ norm_var, scales = "free_y")

nn_analyzing_df %>% 
  filter(Variable != 6 & Variable != 7) %>%
  group_by(Variable, norm_var, norm_error, Method) %>%
  summarize(MSE = mean(MSE)) %>%
  ggplot(aes(x = Variable, y = MSE, fill = Method)) + 
  geom_bar(stat = "identity", position = 'dodge') + 
  facet_wrap(norm_var ~ norm_error, scales = "free_y")

# Summary Tables----

## Bias, Variance and MSE by variable
nn_analyzing_df %>%
  group_by(Variable, Method) %>%
  summarize(Bias = mean(Bias)) %>%
  pivot_wider(names_from = Method, values_from = Bias)

nn_analyzing_df %>%
  group_by(Variable, Method) %>%
  summarize(Variance = mean(Variance)) %>%
  pivot_wider(names_from = Method, values_from = Variance)

nn_analyzing_df %>%
  group_by(Variable, Method) %>%
  summarize(MSE = mean(MSE)) %>%
  pivot_wider(names_from = Method, values_from = MSE)


