## Obtain optimal pH annotations
## Followin the 'mean' approach since the variation isn't that big

library(bugphyzz)
library(dplyr)

ph <- physiologies("optimal ph")[[1]] %>% 
    as_tibble()

ph_summary <- ph %>% 
    select(Parent_name, Attribute_value) %>% 
    group_by(Parent_name) %>% 
    summarise(
        mean = mean(Attribute_value, na.rm = TRUE),
        sd = sd(Attribute_value, na.rm = TRUE),
        min = min(Attribute_value, na.rm = TRUE),
        max = max(Attribute_value, na.rm = TRUE),
        n = n()
    ) %>% 
    ungroup() %>% 
    mutate(
        ph_category = case_when(
            mean >= 0 & mean <= 5.5 ~ 'acidophile',
            mean > 5.5 & mean <= 8 ~ 'neutrophile',
            mean > 8 ~ 'alkaliphile'
        )
    )

ph_annotations <- ph_summary %>% 
    select(Parent_name, ph_category)

# readr::write_tsv(
#     ph_annotations, 
#     file = 'inst/extdata/ph_annotations_bugphyzz.tsv'
# )



