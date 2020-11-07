devtools::load_all()

# setup parallel processing
doParallel::registerDoParallel(parallel::detectCores())
BiocParallel::register(BiocParallel::DoparParam(), default = TRUE)

# raw data from hPMVEC exposed to 21% or 0.5% oxygen for 24 h
# collected on Q Exactive with ZIC-pHILIC separation
# top 10 data-dependent scanning in negative polarity mode
# three 21% and three 0.5% oxygen samples selected as sample dataset

# get full path of data files
files <-
  list.files(path = system.file("extdata", package = "IPO2"),
             pattern = "_1_",
             full.names = TRUE,
             recursive = TRUE)

# create a phenodata data.frame
pd <-
  data.frame(id = stringr::str_extract(basename(files), pattern = "^\\d{2}"),
             group = stringr::str_extract(basename(files), pattern = "(norm|hyp)"),
             stringsAsFactors = FALSE)

# read raw data
raw_data <-
  MSnbase::readMSData(files = files,
                      pdata = new("NAnnotatedDataFrame", pd),
                      msLevel. = 1L,
                      mode = "onDisk")

# perform peak picking with centWave algorithm
cwp <- xcms::CentWaveParam(peakwidth = c(20, 80), noise = 5000)
hpmvec <- xcms::findChromPeaks(raw_data, param = cwp)

usethis::use_data(hpmvec, overwrite = TRUE)
