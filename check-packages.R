requiredPackages <- c(
    "argparse",
    "ggtree",
    "lubridate",
    "scatterpie",
    "tidyr",
    "tidytree",
    "dplyr",
    "config", 
    "DBI"
)

for (p in requiredPackages) {
    if (length(find.package(p, quiet=TRUE)) == 0) {
        print(paste("R package is missing:", p))
        quit(status = 1)
    }
}
