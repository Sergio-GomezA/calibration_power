# update currently used figures in overleaf ####

path1 <- "~/ownCloud-s2441782@datasync.ed.ac.uk/projects/calibration/calibration_power_main_doc/spfig"
path2 <- "~/ownCloud-s2441782@datasync.ed.ac.uk/projects/calibration/calibration_power/fig/batch2025"

# Files to update in overleaf
files1 <- list.files(path1)

# Full paths to matching files in path2
source_files <- file.path(path2, files1)

# Keep only files that actually exist in path2
source_files <- source_files[file.exists(source_files)]

# Destination paths in path1
dest_files <- file.path(path1, basename(source_files))

# Copy and overwrite
file.copy(source_files, dest_files, overwrite = TRUE)
