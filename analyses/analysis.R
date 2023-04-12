library(tidyverse)
library(viridis)
library(scico)
library(ggimage)
library(cmdstanr)
library(bayesplot)
library(tidybayes)
library(bayestestR)

## binning code
source('grid.R')

## set number of cores for parallel computing
options(mc.cores=parallel::detectCores())

fix_encoding <- read_tsv('data/fixreport_encoding2.tsv') %>% mutate(stage='encoding')
fix_simulation <- read_tsv('data/fixreport_simulation2.tsv') %>% mutate(stage='simulation')

fix <- bind_rows(fix_encoding, fix_simulation) %>%
    rename(id=RECORDING_SESSION_LABEL,
           trial=TRIAL_INDEX,
           fix_index=CURRENT_FIX_INDEX,
           fix_x=CURRENT_FIX_X,
           fix_y=CURRENT_FIX_Y,
           fix_start=CURRENT_FIX_START,
           fix_end=CURRENT_FIX_END,
           fix_duration=CURRENT_FIX_DURATION,
           fix_blink=CURRENT_FIX_BLINK_AROUND) %>%
    mutate(id=factor(str_sub(id, end=2)),
           fix_blink=factor(fix_blink),
           decision=factor(decision),
           video=factor(video),
           outcome=factor(outcome),
           prompt=factor(prompt, levels=c('Remember', 'What if?', 'Cause')),
           answer=factor(answer),
           stage=factor(stage),
           ## center fixations around 0
           fix_x=fix_x-1920/2,
           fix_y=-(fix_y-1080/2),
           fix_on_screen=(abs(fix_x) < 1920/2 & abs(fix_y) < 1080/2),
           fix_on_video=(abs(fix_x) < 400 & abs(fix_y) < 300)) %>%
    filter(fix_on_screen & fix_on_video & fix_duration > 75 & fix_index != 1) %>%
    mutate(fix_x=ifelse(decision=='Left', -fix_x, fix_x),   ## reflect fixations so that ball always moves right
           image=ifelse(stage=='encoding', paste0('data/img/', str_to_lower(decision), '_', outcome, '_remember.png'),
                        paste0('data/img/', str_to_lower(decision), '_', outcome, '_', str_to_lower(str_remove_all(prompt, '[\\ \\?]')), '.png')),
           image.reflected=ifelse(stage=='encoding', paste0('data/img/right_', outcome, '_remember.png'),
                                  paste0('data/img/right_', outcome, '_', str_to_lower(str_remove_all(prompt, '[\\ \\?]')), '.png')),
           condition=interaction(stage, prompt, outcome))

img <- fix %>%
    group_by(decision, outcome, stage, prompt) %>%
    summarize(image=image[1],
              image.reflected=image.reflected[1]) %>%
    filter(decision == 'Right')

## pull out behavioral responses
ratings <- fix %>% group_by(id, trial) %>%
    summarize(decision=decision[1], video=video[1], outcome=outcome[1],
              prompt=prompt[1], answer=answer[1])


## make scatterplot
fix %>%
    filter(stage=='encoding') %>%
    ggplot(aes(x=fix_x, y=fix_y, color=id)) +
    geom_image(aes(x=0, y=0, image=image.reflected), size=1,
               by='height', asp=800/600, inherit.aes=FALSE,
               data=img %>% filter(stage=='encoding')) +
    geom_point() +
    coord_fixed(xlim=c(-400, 400), ylim=c(-300, 300), expand=FALSE) +
    facet_grid(outcome ~ prompt) +
    theme_classic() +
    theme(axis.title=element_blank(),
          axis.text=element_blank(),
          axis.ticks=element_blank(),
          axis.line=element_blank())
ggsave('plots/encoding_raw.png', width=8, height=4)

fix %>%
    filter(stage=='simulation') %>%
    ggplot(aes(x=fix_x, y=fix_y, color=id)) +
    geom_image(aes(x=0, y=0, image=image.reflected), size=1,
               by='height', asp=800/600, inherit.aes=FALSE,
               data=img %>% filter(stage=='simulation')) +
    geom_point() +
    coord_fixed(xlim=c(-400, 400), ylim=c(-300, 300), expand=FALSE) +
    facet_grid(outcome ~ prompt) +
    theme_classic() +
    theme(axis.title=element_blank(),
          axis.text=element_blank(),
          axis.ticks=element_blank(),
          axis.line=element_blank())
ggsave('plots/simulation_raw.png', width=8, height=4)


## hexagonal binning
grid_hex <- fix %>% grid.hex(xmin=-400, xmax=400, ymin=-300, ymax=300, xsize=50)
fix_hex <- fix %>% group_by(id, condition, outcome, prompt, stage) %>%
    bin.hex(xmin=-400, xmax=400, ymin=-300, ymax=300, xsize=50)

fix_hex %>%
    filter(stage=='encoding') %>%
    group_by(stage, outcome, prompt, bin_x, bin_y, grid_idx, vi, vx, vy) %>%
    summarize(count=mean(count)) %>%
    unnest(c(vi, vx, vy)) %>%
    ggplot(aes(x=vx, y=vy, group=grid_idx, fill=count, color=count)) +
    geom_polygon(size=.3) +
    scale_color_viridis(option='magma', name='Count', limits=c(0, NA)) +
    scale_fill_viridis(option='magma', name='Count', limits=c(0, NA)) +
    coord_fixed(xlim=c(-400, 400), ylim=c(-300, 300), expand=FALSE) +
    facet_grid(outcome ~ prompt) +
    theme_classic() +
    theme(axis.title=element_blank(),
          axis.text=element_blank(),
          axis.ticks=element_blank(),
          axis.line=element_blank())
ggsave('plots/encoding_hex.png', width=8, height=4)

fix_hex %>%
    filter(stage=='simulation') %>%
    group_by(stage, outcome, prompt, bin_x, bin_y, grid_idx, vi, vx, vy) %>%
    summarize(count=mean(count)) %>%
    unnest(c(vi, vx, vy)) %>%
    ggplot(aes(x=vx, y=vy, group=grid_idx, fill=count, color=count)) +
    geom_polygon(size=.3) +
    scale_color_viridis(option='magma', name='Count', limits=c(0, NA)) +
    scale_fill_viridis(option='magma', name='Count', limits=c(0, NA)) +
    coord_fixed(xlim=c(-400, 400), ylim=c(-300, 300), expand=FALSE) +
    facet_grid(outcome ~ prompt) +
    theme_classic() +
    theme(axis.title=element_blank(),
          axis.text=element_blank(),
          axis.ticks=element_blank(),
          axis.line=element_blank())
ggsave('plots/simulation_hex.png', width=8, height=4)




## model fixations
fix_stan <- list(N=nrow(fix_hex),
                 G=max(fix_hex$grid_idx),
                 P=length(unique(fix_hex$id)),
                 D=2,  ## 2 dimensions (X, Y)
                 C=length(levels(fix_hex$condition)),
                 g=fix_hex$grid_idx,
                 p=as.integer(as.factor(fix_hex$id)),
                 c=as.integer(fix_hex$condition),
                 grid=select(grid_hex, bin_x, bin_y),
                 y=fix_hex$count)
str(fix_stan)

## sample from prior distribution
mod_prior <- cmdstan_model('gp-prior.stan')
if (file.exists('prior_hex_50.rds')) {
    prior <- readRDS('prior_hex_50.rds')
} else {
    prior <- mod_prior$sample(data=fix_stan, iter_warmup=0, iter_sampling=4000, fixed_param=TRUE, output_dir='draws')
    prior$save_object('prior_hex_50.rds')
}

prior$summary(c('prior_a', 'prior_rho', 'prior_rho_tilde',
                'prior_alpha', 'prior_alpha_tilde'))
prior$draws(c('prior_a', 'prior_alpha', 'prior_alpha_tilde')) %>% mcmc_areas(prob=.95)
prior$draws(c('prior_rho', 'prior_rho_tilde')) %>% mcmc_areas(prob=.95)

## sample from posterior distribution
mod <- cmdstan_model('gp.stan')
if (file.exists('gp_hex_50.rds')) {
    fit <- readRDS('gp_hex_50.rds')
} else {
    fit <- mod$sample(data=fix_stan, iter_warmup=1000, iter_sampling=1000, output_dir='draws')
    fit$save_object('gp_hex_50.rds')
}

fit$summary(c('a', 'rho', 'rho_tilde', 'alpha', 'alpha_tilde'))
fit$draws(c('a', 'alpha', 'alpha_tilde')) %>% mcmc_areas(prob=.95)
fit$draws(c('rho', 'rho_tilde')) %>% mcmc_areas(prob=.95)

draws.group <- full_join(spread_draws(prior, prior_a, prior_f[COND, grid_idx]) %>%
                         mutate(.chain=(.iteration-1) %/% 1000 + 1,
                                .iteration=(.iteration-1) %% 1000 + 1),
                         spread_draws(fit, a, f[COND, grid_idx])) %>%
    mutate(prior_lambda=prior_a+prior_f,
           lambda=a+f) %>%
    bind_cols(., grid_hex[.$grid_idx, c('bin_x', 'bin_y', 'vi', 'vx', 'vy')]) %>%
    mutate(condition=levels(fix$condition)[COND]) %>%
    separate(condition, c('stage', 'prompt', 'outcome'), sep='\\.') %>%
    mutate(stage=factor(stage),
           prompt=factor(prompt, levels=c('Remember', 'What if?', 'Cause')),
           outcome=factor(outcome))



draws.group %>%
    filter(stage=='encoding') %>%
    group_by(stage, prompt, outcome, grid_idx, vx, vy) %>%
    median_hdci(lambda) %>%
    unnest(c(vx, vy)) %>%
    ggplot() +
    aes(x=vx, y=vy, group=grid_idx) +
    geom_polygon(aes(fill=exp(lambda), color=exp(lambda)), linewidth=0.3) +
    scale_color_viridis(option='magma', name='Count', limits=c(0, NA)) +
    scale_fill_viridis(option='magma', name='Count', limits=c(0, NA)) +
    coord_fixed(xlim=c(-400, 400), ylim=c(-300, 300), expand=FALSE) +
    facet_grid(outcome ~ prompt) +
    theme_classic() +
    theme(axis.title=element_blank(),
          axis.text=element_blank(),
          axis.ticks=element_blank(),
          axis.line=element_blank())
ggsave('plots/encoding_fit_50.png', width=8, height=4)


draws.group %>%
    filter(stage=='simulation') %>%
    group_by(stage, prompt, outcome, grid_idx, vx, vy) %>%
    median_hdci(lambda) %>%
    unnest(c(vx, vy)) %>%
    ggplot() +
    aes(x=vx, y=vy, group=grid_idx) +
    geom_polygon(aes(fill=exp(lambda), color=exp(lambda)), linewidth=0.3) +
    scale_color_viridis(option='magma', name='Count', limits=c(0, NA)) +
    scale_fill_viridis(option='magma', name='Count', limits=c(0, NA)) +
    coord_fixed(xlim=c(-400, 400), ylim=c(-300, 300), expand=FALSE) +
    facet_grid(outcome ~ prompt) +
    theme_classic() +
    theme(axis.title=element_blank(),
          axis.text=element_blank(),
          axis.ticks=element_blank(),
          axis.line=element_blank())
ggsave('plots/simulation_fit_50.png', width=8, height=4)





## Group-level contrasts
img_contr <- img %>%
    filter(prompt=='Cause') %>%
    select(-prompt) %>%
    expand_grid(prompt=c('What if? - Remember', 'Cause - Remember', 'Cause - What if?')) %>%
    mutate(prompt=factor(prompt, levels=c('What if? - Remember', 'Cause - Remember', 'Cause - What if?')))

contr.group <- draws.group %>%
    group_by(stage, prompt, outcome, grid_idx, bin_x, bin_y, vx, vy) %>%
    compare_levels(lambda, by=prompt) %>%
    left_join(draws.group %>%
              group_by(stage, prompt, outcome, grid_idx, bin_x, bin_y, vx, vy) %>%
              compare_levels(prior_lambda, by=prompt)) %>%
    group_by(stage, prompt, outcome, grid_idx, bin_x, bin_y, vx, vy) %>%
    mutate(log10_BF=bayesfactor_pointnull(lambda, prior_lambda)$log_BF / log(10),
           P=as.numeric(pd_to_p(pd(lambda))),
           prompt=factor(prompt, levels=c('What if? - Remember', 'Cause - Remember', 'Cause - What if?'))) %>%
    group_by(log10_BF, P, .add=TRUE)

contr.group %>%
    filter(log10_BF > 1 & P < .05) %>%    
    median_hdci() %>%
    arrange(stage, prompt, outcome)

## Find significant pixels with the largest effects
clusters.max <- contr.group %>%
    filter(log10_BF > 1 & P < .05) %>%    
    median_hdci() %>%
    filter(lambda > 0) %>%
    group_by(stage, prompt, outcome) %>%
    summarize(i=which.max(lambda), lambda=lambda[i], upper=lambda.upper[i], lower=lambda.lower[i],
              log10_BF=log10_BF[i], P=P[i], grid_idx=grid_idx[i], bin_x=bin_x[i], bin_y=bin_y[i], vx=vx[i], vy=vy[i])
clusters.min <- contr.group %>%
    filter(log10_BF > 1 & P < .05) %>%
    median_hdci() %>%
    filter(lambda < 0) %>%
    group_by(stage, prompt, outcome) %>%
    summarize(i=which.min(lambda), lambda=lambda[i], upper=lambda.upper[i], lower=lambda.lower[i],
              log10_BF=log10_BF[i], P=P[i], grid_idx=grid_idx[i], bin_x=bin_x[i], bin_y=bin_y[i], vx=vx[i], vy=vy[i])
clusters.max %>%
    bind_rows(clusters.min) %>%
    mutate(BF=10^log10_BF) %>%
    relocate(BF, .before=P) %>%
    arrange(stage, prompt, outcome) %>%
    select(-log10_BF, -i, -grid_idx, -vx, -vy) %>%
    write_csv('clusters.csv')


contr.group %>%
    filter(stage=='encoding') %>%
    median_hdci() %>%
    filter(log10_BF > 1 & P < .05) %>%
    unnest(c(vx, vy)) %>%
    ggplot(aes(x=vx, y=vy, group=grid_idx, fill=lambda, color=lambda)) +
    geom_image(aes(x=0, y=0, image=image.reflected), size=1, by='height', asp=800/600, inherit.aes=FALSE,
               data=img_contr %>% filter(stage=='encoding')) +
    geom_polygon(size=0.3, alpha=0.85) +
    scale_fill_scico(palette='vik', midpoint=0, limits=c(-3, 3), breaks=c(-3, -1.5, 0, 1.5, 3), name='Fixation Rate Contrast\n(log scale)',
                     labels=c('-3 (Condition 2 > Condition 1)', '-1.5', '0 (Condition 1 = Condition 2)', '1.5', '3 (Condition 1 > Condition 2)')) +
    scale_color_scico(palette='vik', midpoint=0, limits=c(-3, 3), breaks=c(-3, -1.5, 0, 1.5, 3), name='Fixation Rate Contrast\n(log scale)',
                      labels=c('-3 (Condition 2 > Condition 1)', '-1.5', '0 (Condition 1 = Condition 2)', '1.5', '3 (Condition 1 > Condition 2)')) +
    facet_grid(outcome ~ prompt) +
    coord_fixed(xlim=c(-400, 400), ylim=c(-300, 300), expand=FALSE) +
    theme_classic() +
    theme(axis.title=element_blank(),
          axis.text=element_blank(),
          axis.ticks=element_blank(),
          axis.line=element_blank())
ggsave('plots/encoding_contrast_50_sig.png', width=8, height=4)

contr.group %>%
    filter(stage=='encoding') %>%
    median_hdci() %>%
    unnest(c(vx, vy)) %>%
    ggplot(aes(x=vx, y=vy, group=grid_idx, fill=lambda, color=lambda)) +
    geom_image(aes(x=0, y=0, image=image.reflected), size=1, by='height', asp=800/600, inherit.aes=FALSE,
               data=img_contr %>% filter(stage=='encoding')) +
    geom_polygon(size=0.3, alpha=0.85) +
    scale_fill_scico(palette='vik', midpoint=0, limits=c(-3, 3), breaks=c(-3, -1.5, 0, 1.5, 3), name='Fixation Rate Contrast\n(log scale)',
                     labels=c('-3 (Condition 2 > Condition 1)', '-1.5', '0 (Condition 1 = Condition 2)', '1.5', '3 (Condition 1 > Condition 2)')) +
    scale_color_scico(palette='vik', midpoint=0, limits=c(-3, 3), breaks=c(-3, -1.5, 0, 1.5, 3), name='Fixation Rate Contrast\n(log scale)',
                      labels=c('-3 (Condition 2 > Condition 1)', '-1.5', '0 (Condition 1 = Condition 2)', '1.5', '3 (Condition 1 > Condition 2)')) +
    facet_grid(outcome ~ prompt) +
    coord_fixed(xlim=c(-400, 400), ylim=c(-300, 300), expand=FALSE) +
    theme_classic() +
    theme(axis.title=element_blank(),
          axis.text=element_blank(),
          axis.ticks=element_blank(),
          axis.line=element_blank())
ggsave('plots/encoding_contrast_50.png', width=8, height=4)

contr.group %>%
    filter(stage=='simulation') %>%
    median_hdci() %>%
    filter(log10_BF > 1 & P < .05) %>%
    unnest(c(vx, vy)) %>%
    ggplot(aes(x=vx, y=vy, group=grid_idx, fill=lambda, color=lambda)) +
    geom_image(aes(x=0, y=0, image=image.reflected), size=1, by='height', asp=800/600, inherit.aes=FALSE,
               data=img_contr %>% filter(stage=='simulation')) +
    geom_polygon(size=0.3, alpha=0.85) +
    scale_fill_scico(palette='vik', midpoint=0, limits=c(-3, 3), breaks=c(-3, -1.5, 0, 1.5, 3), name='Fixation Rate Contrast\n(log scale)',
                     labels=c('-3 (Condition 2 > Condition 1)', '-1.5', '0 (Condition 1 = Condition 2)', '1.5', '3 (Condition 1 > Condition 2)')) +
    scale_color_scico(palette='vik', midpoint=0, limits=c(-3, 3), breaks=c(-3, -1.5, 0, 1.5, 3), name='Fixation Rate Contrast\n(log scale)',
                      labels=c('-3 (Condition 2 > Condition 1)', '-1.5', '0 (Condition 1 = Condition 2)', '1.5', '3 (Condition 1 > Condition 2)')) +
    facet_grid(outcome ~ prompt) +
    coord_fixed(xlim=c(-400, 400), ylim=c(-300, 300), expand=FALSE) +
    theme_classic() +
    theme(axis.title=element_blank(),
          axis.text=element_blank(),
          axis.ticks=element_blank(),
          axis.line=element_blank())
ggsave('plots/simulation_contrast_50_sig.png', width=8, height=4)

contr.group %>%
    filter(stage=='simulation') %>%
    median_hdci() %>%
    unnest(c(vx, vy)) %>%
    ggplot(aes(x=vx, y=vy, group=grid_idx, fill=lambda, color=lambda)) +
    geom_image(aes(x=0, y=0, image=image.reflected), size=1, by='height', asp=800/600, inherit.aes=FALSE,
               data=img_contr %>% filter(stage=='simulation')) +
    geom_polygon(size=0.3, alpha=0.85) +
    scale_fill_scico(palette='vik', midpoint=0, limits=c(-3, 3), breaks=c(-3, -1.5, 0, 1.5, 3), name='Fixation Rate Contrast\n(log scale)',
                     labels=c('-3 (Condition 2 > Condition 1)', '-1.5', '0 (Condition 1 = Condition 2)', '1.5', '3 (Condition 1 > Condition 2)')) +
    scale_color_scico(palette='vik', midpoint=0, limits=c(-3, 3), breaks=c(-3, -1.5, 0, 1.5, 3), name='Fixation Rate Contrast\n(log scale)',
                      labels=c('-3 (Condition 2 > Condition 1)', '-1.5', '0 (Condition 1 = Condition 2)', '1.5', '3 (Condition 1 > Condition 2)')) +
    facet_grid(outcome ~ prompt) +
    coord_fixed(xlim=c(-400, 400), ylim=c(-300, 300), expand=FALSE) +
    theme_classic() +
    theme(axis.title=element_blank(),
          axis.text=element_blank(),
          axis.ticks=element_blank(),
          axis.line=element_blank())
ggsave('plots/simulation_contrast_50.png', width=8, height=4)
