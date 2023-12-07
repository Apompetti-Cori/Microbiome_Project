library(tidyverse)
customized_read_tsv <- function(file){
  read_tsv(file, show_col_types = FALSE)
}

suppressMessages(list.files(here("data/metadata/bs_conversion/"),full.names = TRUE) %>% # list all the files
                   lapply(customized_read_tsv) %>% # read them all in with our custom function
                   reduce(bind_rows) -> x)

x <- x %>% 
  select(-H...24, -H_lambda...33) %>% 
  rename(H = H...26, H_lambda = H_lambda...35)

write_tsv(x, file = here("data/metadata/compiled_bs_conversion.tsv"))
