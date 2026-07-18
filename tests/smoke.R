library(RAStools)

result <- get_mean_and_std(c(10, 20, 30), c(0.1, 0.2, 0.3))
stopifnot(result$Nr_blocks == 3)
stopifnot(result$Nr_sites == 60)
stopifnot(isTRUE(all.equal(
    result$f_value,
    sum(c(10, 20, 30) * c(0.1, 0.2, 0.3)) / 60
)))
