#' go to https://www150.statcan.gc.ca/t1/tbl1/en/tv.action?pid=9810040401
#' click add remove data
#' geography==canada
#' age==total
#' highest certificate==select all items
#' major field== all 4 digit
#' noc == 5 digit
#' gender==total
#' download selected for database loading, move to add_to_pond

library(bcgovpond)
ingest_pond()

