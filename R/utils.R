resolve_input_path <- function(path, base_dir = NULL) {
    if (length(path) != 1 || is.na(path) || path == "") {
        return(path)
    }
    if (file.exists(path)) {
        return(path)
    }
    if (!is.null(base_dir)) {
        candidate <- file.path(base_dir, path)
        if (file.exists(candidate)) {
            return(candidate)
        }
    }
    path
}

fill_na_with_ref <- function(tgt, ref) {
    pos <- is.na(tgt)
    tgt[pos] <- if (length(ref) == 1) ref else ref[pos]
    tgt
}
