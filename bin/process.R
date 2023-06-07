#' Pre-process a study with flow cytometry data
#'
#' @return A list of gating set objects.
#' @examples
#' \dontrun{
#' gsl <- process_study()
#' }
process_study <- function() {
  samplesheet <- data.table::fread("samplesheet.csv")
  samplesheet$file_0 <- samplesheet$file_1
  samplesheet$file_1 <- file.path("fcs", samplesheet$file_1)
  samplesheet <- samplesheet[file.exists(samplesheet$file_1)]

  validate_samplesheet(samplesheet)

  files <- summarize_study(samplesheet)

  # create gating set for each panel
  files_by_panel <- split(files, files$panel)
  lapply(seq_along(files_by_panel), function(i) {
    panel <- files_by_panel[[i]]
    gs <- process_panel(panel)
    flowWorkspace::save_gs(gs, paste0("gs", i))
    gs
  })
}

summarize_study <- function(samplesheet) {
  files <- samplesheet

  # summarize files
  catf(sprintf("There are %s fcs files in this dataset", nrow(files)))

  # read headers and summarize panels
  map <- get_channel_map()

  catf("Reading in fcs headers")
  headers <- lapply(files$file_1, function(file) flowCore::read.FCSheader(file, channel_alias = map))

  files$tot <- sapply(headers, function(x) x[[1]]["$TOT"])
  files$par <- sapply(headers, function(x) x[[1]]["$PAR"])

  files$src <- sapply(headers, function(x) x[[1]]["$SRC"])
  files$date <- sapply(headers, function(x) x[[1]]["$DATE"])
  files$btim <- sapply(headers, function(x) x[[1]]["$BTIM"])
  files$etim <- sapply(headers, function(x) x[[1]]["$ETIM"])
  files$cyt <- sapply(headers, function(x) x[[1]]["$CYT"])

  files$creator <- sapply(headers, function(x) x[[1]]["CREATOR"])
  files$tubeName <- sapply(headers, function(x) x[[1]]["TUBE NAME"])
  files$experimentName <- sapply(headers, function(x) x[[1]]["EXPERIMENT NAME"])
  files$settings <- sapply(headers, function(x) x[[1]]["SETTINGS"])
  files$cytnum <- sapply(headers, function(x) x[[1]]["CYTNUM"])
  files$exportUserName <- sapply(headers, function(x) x[[1]]["EXPORT USER NAME"])
  files$exportTime <- sapply(headers, function(x) x[[1]]["EXPORT TIME"])
  files$asf <- sapply(headers, function(x) x[[1]]["FSC ASF"])
  files$plateName <- sapply(headers, function(x) x[[1]]["PLATE NAME"])
  files$plateId <- sapply(headers, function(x) x[[1]]["PLATE ID"])
  files$wellId <- sapply(headers, function(x) x[[1]]["WELL ID"])

  files$cstSetupStatus <- sapply(headers, function(x) x[[1]]["CST SETUP STATUS"])
  files$cstBeadsLotId <- sapply(headers, function(x) x[[1]]["CST BEADS LOT ID"])
  files$cstSetupDate <- sapply(headers, function(x) x[[1]]["CST SETUP DATE"])
  files$cstBaselineDate <- sapply(headers, function(x) x[[1]]["CST BASELINE DATE"])
  files$cytometerConfigName <- sapply(headers, function(x) x[[1]]["CYTOMETER CONFIG NAME"])
  files$cytometerConfigCreateDate <- sapply(headers, function(x) x[[1]]["CYTOMETER CONFIG CREATE DATE"])

  catf("Getting panels for the study")
  panels <- sapply(headers, function(x) {
    header <- x[[1]]
    par <- as.integer(header["$PAR"])
    PNN <- unname(header[paste0("$P", seq_len(par), "N")])
    PNS <- unname(header[paste0("$P", seq_len(par), "S")])

    # standardize channel names
    if (!is.null(map)) {
      for (i in seq_len(nrow(map))) {
        PNN <- gsub(map$channels[i], map$alias[i], PNN)
      }
    }

    PNN <- PNN[!is.na(PNS)]
    PNS <- PNS[!is.na(PNS)]

    if (length(PNS) == 0) {
      ""
    } else {
      paste(gtools::mixedsort(paste0(PNS, " (", PNN, ")")), collapse = "; ")
    }
  })

  files$panel <- panels
  ps <- unique(panels)

  # how many gating sets will be created
  # by panels, sample type, measurement technique, experiment accession
  catf(sprintf("There are %s panel(s)", length(ps)))
  panels_clean <- sapply(strsplit(ps, "; "), function(x) {
    if (length(x) == 0) {
      ""
    } else {
      paste(gsub(" \\(.+\\)$", "", x), collapse = " | ")
    }
  })
  # catf(paste(gtools::mixedsort(panels_clean), collapse = "\n"))

  catf(paste(ps, collapse = "\n"))

  files
}


process_panel <- function(files) {
  panel <- unique(files$panel)
  stopifnot(length(panel) == 1)

  catf(paste(rep("=", times = 80), collapse = ""))
  catf(sprintf("Processing %s fcs files for this panel", nrow(files)))
  catf(paste(strsplit(panel, split = "; ")[[1]], collapse = "\n"))

  # load files
  cs <- create_cytoset(files$file_1)

  # merge metadata
  cs <- merge_metadata(cs, files)

  # create a gating set
  gs <- create_gs(cs)

  # pre-process
  gs <- compensate_gs(gs)
  gs <- transform_gs(gs)

  # gate
  gate_gs(gs)

  gs
}


# processing functions ---------------------------------------------------------
#' @importFrom flowWorkspace load_cytoset_from_fcs cytoset colnames<-
#' @importFrom cytoqc cqc_load_fcs cqc_check cqc_match cqc_match_update cqc_match_remove cqc_fix
create_cytoset <- function(filePath) {
  catf("Reading files and creating a cytoset")
  
  cs <- suppressMessages(cytoqc::cqc_load_fcs(
    filePath,
    num_threads = detect_cores(),
    is_h5 = TRUE
  ))
  channel_check <- cytoqc::cqc_check(cs, "channel")

  # Check if inconsistent
  if (length(unique(channel_check$group_id)) > 1) {
    # Get custom control of channel reference, re-mapping, and fuzzy-matching control from study info in DATA
    #study_info <- DATA[[study]]
    max.distance <- NULL
    channel_ref <- NULL
    map <- get_channel_map()

    channel_match <- custom_match_cytoset(channel_check, max.distance, channel_ref, map)

    # We need to handle possibility of extra channels in reference. By design, cytoqc will not automatically delete these
    # but instead will just throw a warning. For our purposes, if there are extra channels in the reference (e.g. Time, extra scatter channels),
    # we can explicitly delete them to ensure consistency for the resulting cytoset
    missing_channels <- do.call(c, lapply(channel_match$match_result, function(group) {
      group$missing
    }))
    missing_channels <- missing_channels[!missing_channels %in% map$alias]
    if (length(missing_channels) > 0) {
      # Drop those extra channels from the reference
      channel_ref <- channel_match$ref[!channel_match$ref %in% missing_channels]
      # And re-run the match (now the suggested fix will delete them)
      channel_match <- custom_match_cytoset(channel_check, max.distance, channel_ref, map)
    }

    cytoqc::cqc_fix(channel_match)
  }

  # cs with consistent channels
  cs <- flowWorkspace::cytoset(cs)

  # clean scatter channel names
  colnames(cs) <- gsub("^(F|S)S\\d+-(A|W|H)$", "\\1SC-\\2", colnames(cs))

  # clean marker names of non-fluorescence channels
  channels <- c("T0", "T1", "INFO", "FSC-H", "FSC-A", "FSC-W", "SSC-H", "SSC-A", "SSC-W", "TIME")
  to_change <- match(tolower(channels), tolower(colnames(cs)))
  if (any(!is.na(to_change))) {
    channels <- channels[!is.na(to_change)]
    to_change <- to_change[!is.na(to_change)]
    temp_channels <- paste0("TEMP-", channels)
    colnames(cs)[to_change] <- temp_channels
    markernames(cs)[temp_channels] <- NA
    colnames(cs)[to_change] <- channels
  }

  cs
}

# A simple wrapper to handle the HIPCCyto matching logic of:
# 1) Specification of reference channels manually or automatically by most abundant group in check
# 2) Automatic match with optional fuzziness by max.distance
# 3) Manual updates to override automatic match using map
custom_match_cytoset <- function(check_result, max.distance, channel_ref, map) {
  cs <- attr(check_result, "data")
  # 1) If not specified, use the panel with the greatest consensus (most abundant group in check)
  if (is.null(channel_ref)) {
    channel_ref <- colnames(cs[[as.data.frame(check_result)[which.max(check_result$nObject), "object"]]])
  }

  # If not specified, no fuzzy match
  if (is.null(max.distance)) {
    max.distance <- 0.0
  }
  # 2) First try automatic match
  channel_match <- cytoqc::cqc_match(check_result, ref = channel_ref, max.distance = max.distance)

  # 3) Allow manual updating to override automatic match
  if (!is.null(map)) {
    # Remove any existing match
    tryCatch(channel_match <- cytoqc::cqc_match_remove(channel_match, map$channels), error = function(e) {})
    update_ref <- map$alias
    names(update_ref) <- map$channels
    channel_match <- cytoqc::cqc_match_update(channel_match, map = update_ref)
  }
  channel_match
}

#' @importFrom flowWorkspace phenoData phenoData<- cf_keyword_insert
merge_metadata <- function(cs, files) {
  catf("Merging metedata")
  data.table::setkey(files, "file_0")
  files[, name := file_0]
  pd <- Biobase::AnnotatedDataFrame(files[phenoData(cs)$name, ])
  rownames(pd) <- pd$name
  flowWorkspace::phenoData(cs) <- pd

  cs
}

create_gs <- function(cs) {
  catf("Creating a gating set")
  flowWorkspace::GatingSet(cs)
}

#' @importFrom flowWorkspace gs_pop_get_data compensate lapply
#' @importFrom flowCore spillover
compensate_gs <- function(gs, study, debug_dir = NULL) {
  catf("Applying compensation")
  cs <- flowWorkspace::gs_pop_get_data(gs)
  cols <- colnames(gs)
  comp <- lapply(cs, function(x) {
    spills <- flowCore::spillover(x)
    spill <- spills[!sapply(spills, is.null)][[1]] # pick the first non-empty matrix
    keep <- colnames(spill) %in% cols
    spill[keep, keep] # remove extra channels
  })

  gs <- flowWorkspace::compensate(gs, comp)

  gs
}

# transform fluorescence channels using inverse hyperbolic sine transformation
# with cofactor = 150
# https://onlinelibrary.wiley.com/doi/full/10.1002/cyto.a.23030
#' @importFrom flowWorkspace colnames transform transformerList
transform_gs <- function(gs, study, debug_dir = NULL) {
  catf("Applying transformation (inverse hyperbolic sine with cofactor = 150)")
  channels <- colnames2(gs)
  if (length(channels) > 0) {
    cofactor <- 150
    trans_obj <- hipccyto_asinht_trans(cofactor)
    trans_list <- flowWorkspace::transformerList(channels, trans_obj)
    gs <- flowWorkspace::transform(gs, trans_list)
  }

  gs
}

#' @importFrom openCyto gatingTemplate gt_gating gs_add_gating_method
gate_gs <- function(gs) {
  catf("Applying default gating methods")
  # apply_quadrant_gate(gs)
  # apply_singlet_gate(gs, "FSC")
  # apply_singlet_gate(gs, "SSC")
  apply_live_gate(gs)
  apply_nondebris_gate(gs)
  apply_lymphocyte_gate(gs)

  gs
}


# helper functions -------------------------------------------------------------
arcsinh_transform <- function(cofactor, transformationId = "HIPCCytoArcsinh") {
  t <- methods::new("transform", .Data = function(x) asinh(x / cofactor))
  t@transformationId <- transformationId
  t
}

sinh_transform <- function(cofactor, transformationId = "HIPCCytoSinh") {
  t <- methods::new("transform", .Data = function(x) sinh(x) * cofactor)
  t@transformationId <- transformationId
  t
}

hipccyto_asinht_trans <- function(cofactor) {
  trans <- arcsinh_transform(cofactor)
  inv <- sinh_transform(cofactor)
  scales::trans_new("hipccyto_asinht", transform = trans, inverse = inv)
}

get_channel_map <- function() {
  map <- Sys.getenv("CHANNEL_MAP")
  if (!is.null(map)) map <- jsonlite::fromJSON(map)
  map
}

validate_samplesheet <- function(samplesheet) {
  TRUE
}