detect_cores <- function() {
  getOption("mc.cores", parallel::detectCores())
}

catf <- function(..., file = "", sep = " ", fill = TRUE, labels = NULL,
                 append = FALSE) {
  time <- sprintf("[%s]", as.character(Sys.time()))
  cat(time, ..., file = file, sep = sep, fill = fill, labels = labels, append = append)
}

colnames2 <- function(gs) {
  cf <- flowWorkspace::gh_pop_get_data(gs[[1]])
  spillover <- get_spillover(cf)
  colnames(spillover)
}

get_parent <- function(gs) {
  rev(flowWorkspace::gs_get_pop_paths(gs, path = 1))[1]
}

get_spillover <- function(x) {
  spills <- spillover(x)
  spills[!sapply(spills, is.null)][[1]]
}

get_marker_channel <- function(gs) {
  markers <- markernames(gs)
  gsub("/", "_", paste0(markers, "_", names(markers)))
}

get_nodes <- function(gs) {
  gs_get_pop_paths(gs, path = 1)[-1]
}

encode_img <- function(file) {
  paste0("data:", mime::guess_type(file), ";base64,", xfun::base64_encode(file))
}
