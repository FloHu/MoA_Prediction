myhead <- function (x, n = 6L, ncol = 6L, rhs = FALSE) 
{
    # simply a modified version of head.data.frame showing just the first
    # ncol columns
    # argument rhs: if one wants to see the six columns from the right side
    stopifnot(length(n) == 1L)
    stopifnot(length(ncol) == 1L & (ncol > 0))
    n <- if (n < 0L)
        max(nrow(x) + n, 0L)
    else min(n, nrow(x))
    if (!rhs) {
       x[seq_len(n), seq_len(ncol), drop = FALSE]
    } else {
       from <- max(c(ncol(x) - ncol + 1, 1))
       x[seq_len(n), from:(ncol(x))]
    }
}
