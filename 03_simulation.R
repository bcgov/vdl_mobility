library(reticulate)
ot <- import("ot")


#critical to run before any canned transport solver (user is assumed to have ordered things correctly)
stopifnot(identical(names(a), rownames(C)))
stopifnot(identical(names(b), colnames(C)))

P_pot <- ot$sinkhorn(a, b, C, epsilon)




