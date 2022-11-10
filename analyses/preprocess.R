library(tidyverse)

DIR <- 'pilot'

## we'll use this to read in multiple files at once
##d <- read_csv(list.files(paste0('data/', DIR), full.names=TRUE))

## read in just one file for now
d <- read_csv('data/pilot/5sd0muhf6tk4jqetkvwrqaekctfj0dc5.csv') %>%
    filter(id != '' & id != 'simulationbox') %>%
    select(id, stimulus, response, rt) %>%
    mutate(trial=rep(1:(n()/4), each=4)) %>%
    pivot_wider(names_from=id, values_from=c(stimulus, response, rt)) %>%
    select(-response_test, -response_testprompt,
           -rt_test, -rt_testprompt,
           -stimulus_fixationtest) %>%
    // clean up the string columns
    mutate(stimulus_test=str_remove_all(stimulus_test, '\\[|\\]|\\"|(video/)'),
           stimulus_testprompt=str_extract(stimulus_testprompt, '(?<=<p>)(.*)(?=</p>)'),
           stimulus_testquestion=str_remove_all(str_extract(stimulus_testquestion, '(.*)(?=<)'),
                                                '(<p>)|(</p>)|(<.*)'),
           response_fixationtest=str_remove_all(response_fixationtest, 'arrow'),
           response_testquestion=str_remove_all(response_testquestion, 'arrow'))
d %>% select(stimulus_testquestion)

write_csv(d, 'data/processed_data.csv')
