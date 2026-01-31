# temp script to get study IDs

# get files
path <- "data/results/data/binarymat"
files <- list.files(path, recursive = TRUE, full.names = TRUE)
files <- files[-grep("GRAY", files)]

# get filenames
for (file in files) {
    filename <- sub("data/results/data/binarymat/", "", sub("\\.csv", "", file))
    df <- read.csv(file)
    sampleids <- colnames(df)[2:ncol(df)]
    filepath <- paste0(filename, ".txt")
    writeLines(sampleids, con = filepath)
}

file = files[1]
filename <- sub("data/results/data/binarymat/", "", sub("\\.csv", "", file))
df <- read.csv(file)
