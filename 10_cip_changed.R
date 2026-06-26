#evidence that 2016 CIP needs to be mapped... 2 digit changes


library(bcgovpond)
library(tidyverse)
cip <- read_view("CIP-2016-CIP-2021-V1-eng.csv")|>
  mutate(two_2016=str_sub(`CIP 2016 Code / Code CPE 2016`, 1,2),
         two_2021=str_sub(`CIP 2021 V1.0 Code / Code CPE 2021 v1.0`, 1,2),
         two_differ=if_else(two_2016==two_2021, FALSE, TRUE))|>
  select(contains("two"), `CIP 2016 Code / Code CPE 2016`, `CIP 2021 V1.0 Code / Code CPE 2021 v1.0`, contains("title"))|>
  filter(two_differ==TRUE)

