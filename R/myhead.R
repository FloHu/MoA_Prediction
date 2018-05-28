myhead <- function (x, n = 6L, ncol = 6L, ...) 
{
    # simply a modified version of head.data.frame showing just the first
    # ncol columns
    stopifnot(length(n) == 1L)
    stopifnot(length(ncol) == 1L & (ncol > 0))
    n <- if (n < 0L)
        max(nrow(x) + n, 0L)
    else min(n, nrow(x))
    x[seq_len(n), seq_len(ncol), drop = FALSE]
}
