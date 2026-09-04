#!/usr/bin/env Rscript

# Validate DeepRibo's cutoff input and publish its parameter/plot pair only
# after the pinned, checksum-guarded estimator has completed successfully.

usage <- function() {
  cat(paste(
    "Usage: parameter_estimation.R -f DATA_LIST -o PARAMETERS -e ESTIMATOR",
    "[-p S_CURVE_PNG] [-d LEGACY_PLOT_PREFIX] [--receipt COMPLETE]\n"
  ))
}

parse_options <- function(args) {
  values <- list(
    file = NULL,
    out = NULL,
    plot = NULL,
    dest = NULL,
    receipt = NULL,
    engine = NULL
  )
  aliases <- c(
    "-f" = "file", "--file" = "file",
    "-o" = "out", "--out" = "out",
    "-p" = "plot", "--plot" = "plot",
    "-d" = "dest", "--dest" = "dest",
    "--receipt" = "receipt",
    "-e" = "engine", "--engine" = "engine"
  )
  seen <- character()
  index <- 1L
  while (index <= length(args)) {
    argument <- args[[index]]
    if (argument %in% c("-h", "--help")) {
      usage()
      quit(status = 0L)
    }
    key <- unname(aliases[argument])
    if (is.na(key)) {
      stop(paste("unknown option", shQuote(argument)), call. = FALSE)
    }
    if (key %in% seen) {
      stop(paste("option supplied more than once:", argument), call. = FALSE)
    }
    if (index == length(args)) {
      stop(paste("option requires a value:", argument), call. = FALSE)
    }
    values[[key]] <- args[[index + 1L]]
    seen <- c(seen, key)
    index <- index + 2L
  }

  if (is.null(values$file) || is.null(values$out) || is.null(values$engine)) {
    usage()
    stop("input, output, and checksum-guarded estimator paths are required",
         call. = FALSE)
  }
  if (!is.null(values$plot) && !is.null(values$dest)) {
    stop("use either --plot or --dest, not both", call. = FALSE)
  }
  if (is.null(values$plot)) {
    values$plot <- if (is.null(values$dest)) {
      file.path(dirname(values$out), "s_curve.png")
    } else {
      paste0(values$dest, ".png")
    }
  }
  values
}

numeric_column <- function(values, name) {
  numbers <- suppressWarnings(as.numeric(as.character(values)))
  if (length(numbers) != length(values) || any(!is.finite(numbers))) {
    stop(paste("column", shQuote(name), "must contain only finite numbers"),
         call. = FALSE)
  }
  numbers
}

validate_input <- function(path) {
  if (!file.exists(path) || isTRUE(file.info(path)$isdir)) {
    stop(paste("input CSV does not exist:", path), call. = FALSE)
  }
  data <- tryCatch(
    read.csv(path, stringsAsFactors = FALSE, check.names = FALSE),
    error = function(error) {
      stop(paste("cannot read input CSV:", conditionMessage(error)), call. = FALSE)
    }
  )
  required <- c("label", "rpk_elo", "coverage_elo")
  missing <- setdiff(required, names(data))
  if (length(missing) > 0L) {
    stop(paste("input CSV is missing columns:", paste(missing, collapse = ", ")),
         call. = FALSE)
  }

  labels <- numeric_column(data$label, "label")
  rpk <- numeric_column(data$rpk_elo, "rpk_elo")
  coverage <- numeric_column(data$coverage_elo, "coverage_elo")
  if (any(!(labels %in% c(0, 1)))) {
    stop("column 'label' must contain only 0 or 1", call. = FALSE)
  }
  if (any(rpk < 0)) {
    stop("column 'rpk_elo' cannot contain negative occupancy", call. = FALSE)
  }
  if (any(coverage < 0 | coverage > 1)) {
    stop("column 'coverage_elo' must remain within [0, 1]", call. = FALSE)
  }

  usable <- labels == 1 & rpk > 0
  if (sum(usable) < 5L) {
    stop(paste(
      "DeepRibo S-curve estimation requires at least five annotated ORFs",
      "with nonzero elongating-region occupancy"
    ), call. = FALSE)
  }
  if (length(unique(rpk[usable])) < 5L) {
    stop(paste(
      "DeepRibo S-curve estimation requires at least five distinct",
      "positive occupancy values"
    ), call. = FALSE)
  }
  if (length(unique(coverage[usable])) < 2L) {
    stop("DeepRibo S-curve estimation requires variable ORF coverage",
         call. = FALSE)
  }
  invisible(NULL)
}

validate_parameters <- function(parameters) {
  names_required <- c("min_RPKM", "min_coverage")
  if (!is.list(parameters) || !all(names_required %in% names(parameters))) {
    stop("DeepRibo estimator did not return both cutoff values", call. = FALSE)
  }
  values <- vapply(names_required, function(name) {
    value <- suppressWarnings(as.numeric(parameters[[name]]))
    if (length(value) != 1L || !is.finite(value)) {
      stop(paste("DeepRibo estimator returned an invalid", name), call. = FALSE)
    }
    value
  }, numeric(1))
  if (values[["min_RPKM"]] <= 0) {
    stop("DeepRibo estimator returned a non-positive occupancy cutoff",
         call. = FALSE)
  }
  if (values[["min_coverage"]] < 0 || values[["min_coverage"]] > 0.60) {
    stop("DeepRibo estimator returned coverage outside [0, 0.60]",
         call. = FALSE)
  }
  values
}

png_crc32 <- local({
  table <- vapply(0:255, function(value) {
    crc <- as.integer(value)
    for (bit in 1:8) {
      crc <- if (bitwAnd(crc, 1L) != 0L) {
        bitwXor(bitwShiftR(crc, 1L), -306674912L)
      } else {
        bitwShiftR(crc, 1L)
      }
    }
    crc
  }, integer(1))

  function(bytes) {
    crc <- -1L
    for (byte in as.integer(bytes)) {
      index <- bitwAnd(bitwXor(crc, byte), 255L) + 1L
      crc <- bitwXor(bitwShiftR(crc, 8L), table[[index]])
    }
    result <- bitwXor(crc, -1L)
    if (result < 0) as.double(result) + 4294967296 else as.double(result)
  }
})

validate_png <- function(path) {
  link_target <- Sys.readlink(path)
  if ((!is.na(link_target) && nzchar(link_target)) ||
      !isTRUE(file_test("-f", path))) {
    stop("DeepRibo S-curve artifact is not a regular file", call. = FALSE)
  }
  info <- file.info(path)
  if (is.na(info$size) || isTRUE(info$isdir) || info$size < 20L) {
    stop("DeepRibo estimator did not create a non-empty S-curve PNG",
         call. = FALSE)
  }
  data <- readBin(path, what = "raw", n = as.integer(info$size))
  signature <- as.raw(c(0x89, 0x50, 0x4e, 0x47, 0x0d, 0x0a, 0x1a, 0x0a))
  if (!identical(data[seq_along(signature)], signature)) {
    stop("DeepRibo S-curve artifact is not a PNG", call. = FALSE)
  }

  uint32 <- function(bytes) {
    sum(as.double(as.integer(bytes)) * c(16777216, 65536, 256, 1))
  }
  position <- 9L
  first_chunk <- TRUE
  saw_image_data <- FALSE
  saw_end <- FALSE
  while (position <= length(data)) {
    if (position + 11 > length(data)) {
      stop("DeepRibo S-curve artifact is not a PNG", call. = FALSE)
    }
    chunk_length <- uint32(data[position:(position + 3L)])
    chunk_type <- rawToChar(data[(position + 4L):(position + 7L)])
    chunk_end <- position + 11 + chunk_length
    if (!is.finite(chunk_end) || chunk_end > length(data)) {
      stop("DeepRibo S-curve artifact is not a PNG", call. = FALSE)
    }
    checksum_start <- position + 8L + chunk_length
    stored_checksum <- uint32(data[checksum_start:(checksum_start + 3L)])
    checksum_payload <- data[
      (position + 4L):(position + 7L + chunk_length)
    ]
    if (png_crc32(checksum_payload) != stored_checksum) {
      stop("DeepRibo S-curve artifact has an invalid PNG checksum",
           call. = FALSE)
    }
    if (first_chunk) {
      if (chunk_type != "IHDR" || chunk_length != 13) {
        stop("DeepRibo S-curve artifact is not a PNG", call. = FALSE)
      }
      width <- uint32(data[(position + 8L):(position + 11L)])
      height <- uint32(data[(position + 12L):(position + 15L)])
      if (width < 1 || height < 1) {
        stop("DeepRibo S-curve artifact is not a PNG", call. = FALSE)
      }
      first_chunk <- FALSE
    }
    if (chunk_type == "IDAT") {
      saw_image_data <- TRUE
    }
    if (chunk_type == "IEND") {
      if (chunk_length != 0 || chunk_end != length(data)) {
        stop("DeepRibo S-curve artifact is not a PNG", call. = FALSE)
      }
      saw_end <- TRUE
      break
    }
    position <- chunk_end + 1
  }
  if (!saw_image_data || !saw_end) {
    stop("DeepRibo S-curve artifact is not a PNG", call. = FALSE)
  }
  invisible(NULL)
}

path_exists <- function(path) {
  link_target <- Sys.readlink(path)
  file.exists(path) || (!is.na(link_target) && nzchar(link_target))
}

canonical_target_path <- function(path) {
  # normalizePath(..., mustWork = FALSE) returns a missing relative path
  # unchanged, then changes to an absolute path once that output exists.  Build
  # the identity from its existing parent so it is stable across publication.
  file.path(
    normalizePath(dirname(path), mustWork = TRUE),
    basename(path)
  )
}

fsync_paths <- function(files = character(), directories = character()) {
  paths <- c(files, directories)
  if (length(paths) == 0L) {
    return(invisible(NULL))
  }
  python <- Sys.which("python3")
  if (!nzchar(python)) {
    stop("the pinned DeepRibo runtime does not provide python3 for fsync",
         call. = FALSE)
  }
  code <- paste(
    "import os, stat, sys",
    "arguments = iter(sys.argv[1:])",
    "for kind, path in zip(arguments, arguments):",
    "    flags = os.O_RDONLY | getattr(os, 'O_NOFOLLOW', 0)",
    "    if kind == 'directory':",
    "        flags |= getattr(os, 'O_DIRECTORY', 0)",
    "    descriptor = os.open(path, flags)",
    "    try:",
    "        mode = os.fstat(descriptor).st_mode",
    "        if kind == 'directory' and not stat.S_ISDIR(mode):",
    "            raise OSError(path + ' is not a directory')",
    "        if kind == 'file' and not stat.S_ISREG(mode):",
    "            raise OSError(path + ' is not a regular file')",
    "        os.fsync(descriptor)",
    "    finally:",
    "        os.close(descriptor)",
    sep = "\n"
  )
  arguments <- c("-c", shQuote(code))
  for (path in files) {
    arguments <- c(
      arguments,
      "file",
      shQuote(normalizePath(path, mustWork = TRUE))
    )
  }
  for (path in unique(directories)) {
    arguments <- c(
      arguments,
      "directory",
      shQuote(normalizePath(path, mustWork = TRUE))
    )
  }
  status <- system2(python, arguments, stdout = "", stderr = "")
  if (status != 0L) {
    stop("cannot synchronize DeepRibo publication state", call. = FALSE)
  }
  invisible(NULL)
}

validate_existing_target <- function(path) {
  if (!path_exists(path)) {
    return(invisible(NULL))
  }
  link_target <- Sys.readlink(path)
  if ((!is.na(link_target) && nzchar(link_target)) ||
      !isTRUE(file_test("-f", path))) {
    stop(paste("output path is not a regular file:", path), call. = FALSE)
  }
  invisible(NULL)
}

remove_file <- function(path, label) {
  if (!path_exists(path)) {
    return(invisible(NULL))
  }
  if (isTRUE(file.info(path)$isdir)) {
    stop(paste(label, "is a directory:", path), call. = FALSE)
  }
  status <- unlink(path, recursive = FALSE, force = TRUE)
  if (status != 0L || path_exists(path)) {
    stop(paste("cannot remove", label, path), call. = FALSE)
  }
  invisible(NULL)
}

remove_transaction <- function(path) {
  if (!path_exists(path)) {
    return(invisible(NULL))
  }
  link_target <- Sys.readlink(path)
  if ((!is.na(link_target) && nzchar(link_target)) ||
      !isTRUE(file.info(path)$isdir)) {
    stop(paste("publication transaction is not a real directory:", path),
         call. = FALSE)
  }
  status <- unlink(path, recursive = TRUE, force = TRUE)
  if (status != 0L || path_exists(path)) {
    stop(paste("cannot remove publication transaction:", path),
         call. = FALSE)
  }
  invisible(NULL)
}

mark_transaction_cleanup_only <- function(path) {
  ready <- file.path(path, "ready")
  if (!transaction_marker_exists(ready, "ready marker")) {
    return(invisible(NULL))
  }
  remove_file(ready, "publication transaction ready marker")
  # The absence of `ready` means that all output restoration/commit work is
  # complete and only disposable transaction metadata remains.  Persist that
  # transition before recursively removing the directory: a later invocation
  # can then safely discard even a partially deleted marker set.
  fsync_paths(directories = path)
  invisible(NULL)
}

transaction_paths <- function(targets) {
  # Two output paths cannot be replaced with one POSIX rename. Persist the
  # publication phase beside the parameter file so a later invocation can
  # distinguish an interrupted replacement from a completed pair.
  parent <- dirname(canonical_target_path(targets[[1L]]))
  prefix <- paste0(".", basename(targets[[1L]]), ".pair-transaction")
  list(
    staging = file.path(parent, paste0(prefix, ".staging")),
    committed = file.path(parent, paste0(prefix, ".committed"))
  )
}

transaction_member <- function(transaction, kind, index) {
  file.path(transaction, paste0(kind, "-", index))
}

transaction_backup <- function(target, index) {
  file.path(
    dirname(target),
    paste0(".", basename(target), ".pair-transaction-old-", index)
  )
}

create_transaction_marker <- function(path, content = NULL) {
  if (!isTRUE(file.create(path))) {
    stop(paste("cannot create publication transaction marker:", path),
         call. = FALSE)
  }
  if (!is.null(content)) {
    writeLines(content, con = path, useBytes = TRUE)
  }
  invisible(NULL)
}

receipt_temporary <- function(path) {
  paste0(path, ".tmp")
}

invalidate_receipt <- function(path) {
  temporary <- receipt_temporary(path)
  validate_existing_target(path)
  validate_existing_target(temporary)
  remove_file(temporary, "temporary cutoff receipt")
  remove_file(path, "cutoff receipt")
  # Absence is the durable signal that no consumer may use the pair while a
  # runner is reconciling or publishing it.
  fsync_paths(directories = dirname(path))
  invisible(NULL)
}

publish_receipt <- function(path, targets) {
  temporary <- receipt_temporary(path)
  if (path_exists(path) || path_exists(temporary)) {
    stop("cutoff publication receipt already exists", call. = FALSE)
  }
  create_transaction_marker(
    temporary,
    c(
      "HRIBO DeepRibo cutoff publication v1",
      unname(vapply(targets, canonical_target_path, character(1)))
    )
  )
  fsync_paths(files = temporary, directories = dirname(temporary))
  if (!file.rename(temporary, path)) {
    stop("cannot publish cutoff completion receipt", call. = FALSE)
  }
  fsync_paths(files = path, directories = dirname(path))
  invisible(NULL)
}

transaction_marker_exists <- function(path, label) {
  if (!path_exists(path)) {
    return(FALSE)
  }
  link_target <- Sys.readlink(path)
  if ((!is.na(link_target) && nzchar(link_target)) ||
      !isTRUE(file_test("-f", path))) {
    stop(paste("publication transaction", label,
               "is not a regular file:", path), call. = FALSE)
  }
  TRUE
}

inspect_transaction <- function(transaction_path, targets) {
  link_target <- Sys.readlink(transaction_path)
  if ((!is.na(link_target) && nzchar(link_target)) ||
      !isTRUE(file.info(transaction_path)$isdir)) {
    stop(paste("publication transaction is not a real directory:",
               transaction_path), call. = FALSE)
  }

  ready <- transaction_marker_exists(
    file.path(transaction_path, "ready"), "ready marker"
  )
  if (ready) {
    identity_paths <- vapply(
      seq_along(targets),
      function(index) transaction_member(transaction_path, "target", index),
      character(1)
    )
    identities <- vapply(
      seq_along(identity_paths),
      function(index) transaction_marker_exists(
        identity_paths[[index]], paste("target marker", index)
      ),
      logical(1)
    )
    if (!all(identities)) {
      stop("publication transaction has incomplete target identity",
           call. = FALSE)
    }
    for (index in seq_along(targets)) {
      stored <- readLines(identity_paths[[index]], warn = FALSE)
      expected <- canonical_target_path(targets[[index]])
      if (length(stored) != 1L || !identical(stored[[1L]], expected)) {
        stop(paste("publication transaction target mismatch for output", index),
             call. = FALSE)
      }
    }
  }
  states <- vector("list", length(targets))
  for (index in seq_along(targets)) {
    present <- transaction_marker_exists(
      transaction_member(transaction_path, "present", index),
      paste("present marker", index)
    )
    absent <- transaction_marker_exists(
      transaction_member(transaction_path, "absent", index),
      paste("absent marker", index)
    )
    backup <- transaction_backup(targets[[index]], index)
    backup_exists <- path_exists(backup)
    if (present && absent) {
      stop(paste("publication transaction has conflicting state for output",
                 index), call. = FALSE)
    }
    if (!ready && backup_exists) {
      stop("publication transaction was modified before it became ready",
           call. = FALSE)
    }
    if (ready && !xor(present, absent)) {
      stop(paste("publication transaction has incomplete state for output",
                 index), call. = FALSE)
    }
    if (absent && backup_exists) {
      stop(paste("publication transaction has a backup for absent output",
                 index), call. = FALSE)
    }
    backup_link <- Sys.readlink(backup)
    if (backup_exists &&
        ((!is.na(backup_link) && nzchar(backup_link)) ||
         !isTRUE(file_test("-f", backup)))) {
      stop(paste("publication transaction backup is not a regular file:", backup),
           call. = FALSE)
    }
    if (ready && present && !backup_exists && !path_exists(targets[[index]])) {
      stop(paste("publication transaction lost existing output", index),
           call. = FALSE)
    }
    states[[index]] <- list(
      present = present,
      absent = absent,
      backup = backup,
      backup_exists = backup_exists
    )
  }
  list(ready = ready, states = states)
}

rollback_staging_transaction <- function(staging, targets) {
  transaction <- inspect_transaction(staging, targets)
  if (!transaction$ready) {
    remove_transaction(staging)
    fsync_paths(directories = dirname(staging))
    return(invisible(NULL))
  }

  for (index in rev(seq_along(targets))) {
    state <- transaction$states[[index]]
    target <- targets[[index]]
    if (state$backup_exists) {
      remove_file(target, "partially published output")
      if (!file.rename(state$backup, target)) {
        stop(paste("cannot restore interrupted output:", target),
             call. = FALSE)
      }
    } else if (state$absent) {
      remove_file(target, "partially published output")
    }
  }
  existing_targets <- targets[vapply(targets, path_exists, logical(1))]
  fsync_paths(
    files = existing_targets,
    directories = unique(dirname(targets))
  )
  mark_transaction_cleanup_only(staging)
  remove_transaction(staging)
  fsync_paths(directories = dirname(staging))
  invisible(NULL)
}

validate_parameter_artifact <- function(path) {
  link_target <- Sys.readlink(path)
  if (!file.exists(path) || (!is.na(link_target) && nzchar(link_target)) ||
      !isTRUE(file_test("-f", path))) {
    stop("committed publication is missing its parameter output",
         call. = FALSE)
  }
  lines <- tryCatch(
    readLines(path, warn = FALSE),
    error = function(error) {
      stop(paste("cannot read committed parameter output:",
                 conditionMessage(error)), call. = FALSE)
    }
  )
  fields <- if (length(lines) == 1L) {
    strsplit(lines[[1L]], ",", fixed = TRUE)[[1L]]
  } else {
    character()
  }
  values <- suppressWarnings(as.numeric(fields))
  if (length(fields) != 2L || any(!nzchar(trimws(fields))) ||
      length(values) != 2L || any(!is.finite(values)) || values[[1L]] <= 0 ||
      values[[2L]] < 0 || values[[2L]] > 0.60) {
    stop("committed publication has an invalid parameter output",
         call. = FALSE)
  }
  invisible(NULL)
}

legacy_backups <- function(target) {
  directory <- dirname(target)
  prefix <- paste0(".", basename(target), ".previous.")
  entries <- list.files(
    directory,
    all.files = TRUE,
    full.names = TRUE,
    no.. = TRUE
  )
  sort(entries[startsWith(basename(entries), prefix)])
}

validate_backup_file <- function(path, label) {
  link_target <- Sys.readlink(path)
  if ((!is.na(link_target) && nzchar(link_target)) ||
      !isTRUE(file_test("-f", path))) {
    stop(paste(label, "is not a regular file:", path), call. = FALSE)
  }
  invisible(NULL)
}

legacy_recovery_paths <- function(targets) {
  parent <- dirname(canonical_target_path(targets[[1L]]))
  marker <- file.path(
    parent,
    paste0(".", basename(targets[[1L]]), ".pair-legacy-restore")
  )
  list(marker = marker, temporary = paste0(marker, ".tmp"))
}

inspect_legacy_recovery <- function(targets) {
  paths <- legacy_recovery_paths(targets)
  active <- transaction_marker_exists(
    paths$marker, "legacy-recovery marker"
  )
  temporary <- transaction_marker_exists(
    paths$temporary, "temporary legacy-recovery marker"
  )
  if (active && temporary) {
    stop("both legacy-recovery marker states exist; refusing to guess",
         call. = FALSE)
  }
  if (temporary) {
    # Output mutation starts only after the completed marker has been renamed
    # into place and synchronized, so an orphan temporary marker is disposable.
    remove_file(paths$temporary, "temporary legacy-recovery marker")
    fsync_paths(directories = dirname(paths$temporary))
  }
  if (active) {
    stored <- readLines(paths$marker, warn = FALSE)
    expected <- unname(vapply(
      targets, canonical_target_path, character(1)
    ))
    if (!identical(stored, expected)) {
      stop("legacy-recovery marker target mismatch", call. = FALSE)
    }
  }
  list(paths = paths, active = active)
}

start_legacy_recovery <- function(targets, paths) {
  if (path_exists(paths$marker) || path_exists(paths$temporary)) {
    stop("legacy-recovery marker already exists", call. = FALSE)
  }
  create_transaction_marker(
    paths$temporary,
    unname(vapply(targets, canonical_target_path, character(1)))
  )
  fsync_paths(
    files = paths$temporary,
    directories = dirname(paths$temporary)
  )
  if (!file.rename(paths$temporary, paths$marker)) {
    stop("cannot publish legacy-recovery marker", call. = FALSE)
  }
  fsync_paths(directories = dirname(paths$marker))
  invisible(NULL)
}

finish_legacy_recovery <- function(paths) {
  remove_file(paths$marker, "legacy-recovery marker")
  fsync_paths(directories = dirname(paths$marker))
  invisible(NULL)
}

restore_legacy_snapshot <- function(targets, backups, recovery_paths) {
  # Keep both source backups until the complete restored pair has been
  # validated and synchronized.  If the process dies (or the second copy
  # fails), the next invocation can repeat every copy from the intact pair.
  for (index in seq_along(targets)) {
    target <- targets[[index]]
    source <- backups[[index]][[1L]]
    remove_file(target, "partially published legacy output")
    copied <- suppressWarnings(file.copy(
      source,
      target,
      overwrite = FALSE,
      copy.mode = TRUE,
      copy.date = TRUE
    ))
    if (!isTRUE(copied)) {
      stop(paste("cannot restore legacy publication output:", target),
           call. = FALSE)
    }
    validate_existing_target(target)
  }

  validate_parameter_artifact(targets[[1L]])
  validate_png(targets[[2L]])
  fsync_paths(files = targets, directories = unique(dirname(targets)))
  # Once this deletion is durable, a surviving backup means only that cleanup
  # was interrupted; before it, both backups remain authoritative rollback
  # sources even if a partly restored pair happens to validate individually.
  finish_legacy_recovery(recovery_paths)
  for (paths in backups) {
    remove_file(paths[[1L]], "legacy publication backup")
  }
  fsync_paths(directories = unique(dirname(targets)))
  invisible(NULL)
}

reconcile_legacy_backups <- function(targets) {
  recovery <- inspect_legacy_recovery(targets)
  backups <- lapply(targets, legacy_backups)
  if (any(lengths(backups) > 1L)) {
    stop("multiple legacy publication backups exist; refusing to guess",
         call. = FALSE)
  }
  if (!any(lengths(backups) == 1L) && !recovery$active) {
    return(invisible(NULL))
  }
  for (paths in backups) {
    if (length(paths) == 1L) {
      validate_backup_file(paths[[1L]], "legacy publication backup")
    }
  }

  if (recovery$active) {
    if (!all(lengths(backups) == 1L)) {
      stop("active legacy restoration has an incomplete snapshot",
           call. = FALSE)
    }
    restore_legacy_snapshot(targets, backups, recovery$paths)
    return(invisible(NULL))
  }

  current <- vapply(targets, path_exists, logical(1))
  if (any(!current)) {
    if (!all(lengths(backups) == 1L)) {
      stop("legacy publication snapshot is incomplete; refusing to guess",
           call. = FALSE)
    }
    start_legacy_recovery(targets, recovery$paths)
    restore_legacy_snapshot(targets, backups, recovery$paths)
    return(invisible(NULL))
  }

  validation_error <- tryCatch(
    {
      validate_parameter_artifact(targets[[1L]])
      validate_png(targets[[2L]])
      NULL
    },
    error = function(error) conditionMessage(error)
  )
  if (!is.null(validation_error)) {
    if (!all(lengths(backups) == 1L)) {
      stop(paste(
        "legacy publication is invalid and its previous pair is incomplete:",
        validation_error
      ), call. = FALSE)
    }
    start_legacy_recovery(targets, recovery$paths)
    restore_legacy_snapshot(targets, backups, recovery$paths)
    return(invisible(NULL))
  } else {
    for (paths in backups) {
      if (length(paths) == 1L) {
        remove_file(paths[[1L]], "legacy publication backup")
      }
    }
  }
  fsync_paths(files = targets, directories = unique(dirname(targets)))
  invisible(NULL)
}

reconcile_publication <- function(
    targets, transactions, defer_committed_cleanup_failure = FALSE) {
  has_staging <- path_exists(transactions$staging)
  has_committed <- path_exists(transactions$committed)
  if (has_staging && has_committed) {
    stop("both staging and committed publication transactions exist",
         call. = FALSE)
  }
  if (!has_staging && !has_committed) {
    reconcile_legacy_backups(targets)
    deterministic_backups <- vapply(
      seq_along(targets),
      function(index) path_exists(transaction_backup(targets[[index]], index)),
      logical(1)
    )
    if (any(deterministic_backups)) {
      stop("orphaned publication backup exists without transaction state",
           call. = FALSE)
    }
  }
  if (has_staging) {
    rollback_staging_transaction(transactions$staging, targets)
  }
  if (has_committed) {
    transaction <- inspect_transaction(
      transactions$committed,
      targets
    )
    if (!transaction$ready) {
      # A durably cleared ready marker is the cleanup-only state.  Marker files
      # may be missing because a prior recursive deletion was interrupted.
      remove_transaction(transactions$committed)
      fsync_paths(directories = dirname(transactions$committed))
      return(invisible(NULL))
    }
    validation_error <- tryCatch(
      {
        validate_parameter_artifact(targets[[1L]])
        validate_png(targets[[2L]])
        NULL
      },
      error = function(error) conditionMessage(error)
    )
    if (!is.null(validation_error)) {
      complete_snapshot <- all(vapply(
        transaction$states,
        function(state) state$absent || state$backup_exists,
        logical(1)
      ))
      if (complete_snapshot) {
        rollback_staging_transaction(transactions$committed, targets)
        return(invisible(NULL))
      }
      stop(paste(
        "committed publication is invalid and its previous pair is incomplete:",
        validation_error
      ), call. = FALSE)
    }
    cleanup_error <- tryCatch(
      {
        for (state in transaction$states) {
          if (state$backup_exists) {
            remove_file(state$backup, "committed publication backup")
          }
        }
        fsync_paths(directories = unique(dirname(targets)))
        mark_transaction_cleanup_only(transactions$committed)
        remove_transaction(transactions$committed)
        fsync_paths(directories = dirname(transactions$committed))
        NULL
      },
      error = function(error) conditionMessage(error)
    )
    if (!is.null(cleanup_error)) {
      if (!defer_committed_cleanup_failure) {
        stop(cleanup_error, call. = FALSE)
      }
      cat(
        "parameter_estimation: warning: committed outputs are valid; ",
        "transaction cleanup is deferred: ", cleanup_error, "\n",
        sep = "", file = stderr()
      )
    }
  }
  invisible(NULL)
}

publish_pair <- function(temporary, targets) {
  transactions <- transaction_paths(targets)
  reconcile_publication(targets, transactions)
  for (target in targets) {
    validate_existing_target(target)
  }

  tryCatch({
    if (!dir.create(transactions$staging, recursive = FALSE)) {
      stop(paste("cannot create publication transaction:",
                 transactions$staging), call. = FALSE)
    }
    for (index in seq_along(targets)) {
      target <- targets[[index]]
      backup <- transaction_backup(target, index)
      if (path_exists(backup)) {
        stop(paste("orphaned publication backup already exists:", backup),
             call. = FALSE)
      }
      create_transaction_marker(
        transaction_member(transactions$staging, "target", index),
        canonical_target_path(target)
      )
      if (path_exists(target)) {
        create_transaction_marker(
          transaction_member(transactions$staging, "present", index)
        )
      } else {
        create_transaction_marker(
          transaction_member(transactions$staging, "absent", index)
        )
      }
    }
    transaction_state_files <- list.files(
      transactions$staging,
      all.files = TRUE,
      full.names = TRUE,
      no.. = TRUE
    )
    fsync_paths(
      files = transaction_state_files,
      directories = c(
        transactions$staging,
        dirname(transactions$staging)
      )
    )
    # `ready` authorizes output mutation. Persist the complete target/state
    # description first, then publish and synchronize this barrier separately.
    ready <- file.path(transactions$staging, "ready")
    create_transaction_marker(ready)
    fsync_paths(files = ready, directories = transactions$staging)

    for (index in seq_along(targets)) {
      target <- targets[[index]]
      present <- transaction_member(transactions$staging, "present", index)
      if (path_exists(present) &&
          !file.rename(
            target,
            transaction_backup(target, index)
          )) {
          stop(paste("cannot stage existing output for replacement:", target),
               call. = FALSE)
      }
    }
    fsync_paths(directories = unique(dirname(targets)))
    for (index in seq_along(targets)) {
      if (!file.rename(temporary[[index]], targets[[index]])) {
        stop(paste("cannot publish output:", targets[[index]]), call. = FALSE)
      }
    }
    fsync_paths(files = targets, directories = unique(dirname(targets)))
    if (!file.rename(transactions$staging, transactions$committed)) {
      stop("cannot commit parameter and S-curve publication", call. = FALSE)
    }
    fsync_paths(directories = dirname(transactions$staging))
  }, error = function(error) {
    rollback_error <- tryCatch(
      {
        if (path_exists(transactions$staging)) {
          rollback_staging_transaction(transactions$staging, targets)
        }
        NULL
      },
      error = function(rollback_condition) conditionMessage(rollback_condition)
    )
    if (!is.null(rollback_error)) {
      stop(paste(conditionMessage(error),
                 "rollback also failed:", rollback_error), call. = FALSE)
    }
    stop(conditionMessage(error), call. = FALSE)
  })

  reconcile_publication(
    targets, transactions, defer_committed_cleanup_failure = TRUE
  )
  invisible(NULL)
}

main <- function() {
  options <- parse_options(commandArgs(trailingOnly = TRUE))
  output_directories <- unique(c(
    dirname(options$out),
    dirname(options$plot),
    if (!is.null(options$receipt)) dirname(options$receipt) else character()
  ))
  for (directory in output_directories) {
    if (!dir.exists(directory) && !dir.create(directory, recursive = TRUE)) {
      stop(paste("cannot create output directory:", directory), call. = FALSE)
    }
  }
  output_path <- canonical_target_path(options$out)
  plot_path <- canonical_target_path(options$plot)
  receipt_path <- if (!is.null(options$receipt)) {
    canonical_target_path(options$receipt)
  } else {
    NULL
  }
  if (identical(output_path, plot_path)) {
    stop("parameter and S-curve outputs must use different paths", call. = FALSE)
  }
  if (!is.null(receipt_path) &&
      receipt_path %in% c(output_path, plot_path)) {
    stop("cutoff receipt must use a different path from its artifacts",
         call. = FALSE)
  }
  protected_paths <- normalizePath(
    c(options$file, options$engine),
    mustWork = FALSE
  )
  if (any(c(output_path, plot_path, receipt_path) %in% protected_paths)) {
    stop("cutoff outputs must not overwrite an input or estimator",
         call. = FALSE)
  }
  targets <- c(output_path, plot_path)
  if (!is.null(receipt_path)) {
    invalidate_receipt(receipt_path)
  }
  transactions <- transaction_paths(targets)
  reconcile_publication(targets, transactions)
  for (target in targets) {
    validate_existing_target(target)
  }
  if (!file.exists(options$file) || isTRUE(file.info(options$file)$isdir)) {
    stop(paste("input CSV does not exist:", options$file), call. = FALSE)
  }
  if (!file.exists(options$engine) || isTRUE(file.info(options$engine)$isdir)) {
    stop(paste("S-curve estimator does not exist:", options$engine), call. = FALSE)
  }
  validate_input(options$file)

  parameter_temporary <- tempfile(
    pattern = paste0(".", basename(output_path), "."),
    tmpdir = dirname(output_path)
  )
  plot_prefix <- tempfile(
    pattern = paste0(".", basename(plot_path), "."),
    tmpdir = dirname(plot_path)
  )
  plot_temporary <- paste0(plot_prefix, ".png")
  on.exit(unlink(c(parameter_temporary, plot_temporary), force = TRUE), add = TRUE)

  estimator_environment <- new.env(parent = globalenv())
  tryCatch(
    sys.source(options$engine, envir = estimator_environment),
    error = function(error) {
      stop(paste("cannot load S-curve estimator:", conditionMessage(error)),
           call. = FALSE)
    }
  )
  if (!exists("get_cutoff_values", envir = estimator_environment,
              mode = "function", inherits = FALSE)) {
    stop("S-curve source does not define get_cutoff_values", call. = FALSE)
  }

  parameters <- tryCatch(
    estimator_environment$get_cutoff_values(
      path = options$file,
      dest = plot_prefix
    ),
    error = function(error) {
      stop(paste("DeepRibo S-curve estimation failed:",
                 conditionMessage(error)), call. = FALSE)
    }
  )
  values <- validate_parameters(parameters)
  validate_png(plot_temporary)

  output <- paste(
    sprintf("%.17g", values[["min_RPKM"]]),
    sprintf("%.17g", values[["min_coverage"]]),
    sep = ","
  )
  writeLines(output, con = parameter_temporary, useBytes = TRUE)
  fsync_paths(
    files = c(parameter_temporary, plot_temporary),
    directories = unique(dirname(c(parameter_temporary, plot_temporary)))
  )
  publish_pair(
    c(parameter_temporary, plot_temporary),
    targets
  )
  if (!is.null(receipt_path)) {
    publish_receipt(receipt_path, targets)
  }
}

tryCatch(
  main(),
  error = function(error) {
    cat("parameter_estimation: error:", conditionMessage(error), "\n",
        file = stderr())
    quit(status = 1L)
  }
)
