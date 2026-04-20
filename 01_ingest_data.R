#' go to https://www150.statcan.gc.ca/t1/tbl1/en/tv.action?pid=9810040401
#' click add remove data
#' geography==canada
#' age==total
#' highest certificate==select all items
#' major field== all 4 digit
#' noc == 5 digit
#' gender==total
#' download selected for database loading, move to add_to_pond
#'
#' ONET: go to https://www.onetcenter.org/dl_files/database/db_30_2_excel.zip
#' unzip file, move skills, abilities, knowledge, work activities to add_to_pond
#' prepend each file with date and "_" e.g. 2026-02-16_abilities.xlsx  it is crucial there is no "_" in the date, and
#' the date and the meaning of the file are separated by "_"  (bcgovpond)
#'
#' then run following:

library(bcgovpond)
ingest_pond()

