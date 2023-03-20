library(tidyverse)
library(viridis)
library(scico)
library(ggimage)


fix_encoding <- read_tsv('data/fixreport_encoding.tsv') %>% mutate(stage='encoding')
fix_simulation <- read_tsv('data/fixreport_simulation.tsv') %>% mutate(stage='simulation')

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
           image=paste0('data/img/LAMI_kids_', str_to_lower(decision), '_', outcome, '.png'),
           image.reflected=paste0('data/img/LAMI_kids_right_', outcome, '.png'))

img <- fix %>%
    group_by(decision, outcome, stage) %>%
    summarize(image=image[1],
              image.reflected=image.reflected[1]) %>%
    filter(decision == 'Right')



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
source('grid.R')

fix_hex <- fix %>% group_by(id, outcome, prompt, stage) %>%
    bin.hex(xmin=-400, xmax=400, ymin=-300, ymax=300, xsize=60)

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


## rectangular binning
fix_rect <- fix %>% group_by(id, outcome, prompt, stage) %>%
    bin.rect(xmin=-200, xmax=200, xsize=400,
             ymin=-150, ymax=150, ysize=300)

fix_rect %>%
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
ggsave('plots/encoding_rect.png', width=8, height=4)

fix_rect %>%
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
ggsave('plots/simulation_rect.png', width=8, height=4)




## pull out behavioral responses
ratings <- fix %>% group_by(id, trial) %>%
    summarize(decision=decision[1], video=video[1], outcome=outcome[1],
              prompt=prompt[1], answer=answer[1])
