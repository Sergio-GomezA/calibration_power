# update currently used figures in overleaf ####

path1 <- "~/ownCloud-s2441782@datasync.ed.ac.uk/projects/calibration/calibration_power_main_doc/spfig"
path2 <- "~/ownCloud-s2441782@datasync.ed.ac.uk/projects/calibration/calibration_power/fig/batchY25d150_v5"
path_alt <- c(
  "/home/s2441782/Documents/elexon/caloutput/batchY25d150_v5/fig/oos",
  "/home/s2441782/Documents/elexon/caloutput/batchY25d150_v5/fig/fit",
  "/home/s2441782/Documents/elexon/caloutput/batchY25d150_v5/fig/oos/by_model"
)

# Files to update in overleaf
files1 <- list.files(path1)

# Full paths to matching files in path2
source_files <- file.path(path2, files1)

# Keep only files that actually exist in path2
source_files <- source_files[file.exists(source_files)]

# Search all alternative folders
for (p in path_alt) {
  alt_files <- file.path(p, files1)
  alt_files <- alt_files[file.exists(alt_files)]

  source_files <- c(source_files, alt_files)
}

# Copy all matching files to path1
dest_files <- file.path(path1, basename(source_files))

# Copy and overwrite
file.copy(source_files, dest_files, overwrite = TRUE)
