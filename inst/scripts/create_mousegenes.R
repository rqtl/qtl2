# create sqlite database with mouse genes, from MGI
#
# source (dated 2017-11-03):
#      http://www.informatics.jax.org/downloads/mgigff/
#          MGI.20171103.gff3.gz
#          MGI_GFF_Spec.docx (annotations, e.g. column names)

library(RSQLite)

### download files
site <- "http://www.informatics.jax.org/downloads/mgigff"
file <- "MGI.20171103.gff3.gz"
date_source <- ymd(strsplit(file, "\\.")[[1]][2])
genome_build <- "GRCm38/mm10"
url <- paste0(site, "/", file)
if(!file.exists(file))
    download.file(url, file)

# load gene table
file <- sub(".gz$", "", file)
tmpfile <- tempfile()
remove_tmpfile <- FALSE
if(!file.exists(file)) { # need to unzip
    system(paste0("gunzip -c ", file, ".gz > ", tmpfile))
    remove_tmpfile <- TRUE
    file <- tmpfile
}

# fields in the file
fields <- c("chr", "source", "type", "start", "stop", "score", "strand",
            "phase", "attributes")
# within attributes field:
attrib <- c("ID", "Name", "Parent", "Dbxref", "mgiName", "bioType", "Alias")

# read data
tab <- read.table(file, sep="\t", header=FALSE, comment.char="#",
                  na.strings=".", stringsAsFactors=FALSE,
                  quote="", fill=FALSE,
                  colClasses=rep(rep(c("character", "numeric"), 3), c(3,3,1,1,1,0)))
colnames(tab) <- fields

# remove temporary file
if(remove_tmpfile) {
    unlink(tmpfile)
}

# remove things not on chr 1-19, "X", "Y", "MT"
tab <- tab[tab$chr %in% c(1:19, "X", "Y", "MT"),]

# split 9th column at ';' then split at '=' and use first bit as key and second bit as value
tab9_spl <- strsplit(tab[,9], ";")
tab9_list <- lapply(tab9_spl, function(a) {
    spl <- strsplit(a, "=", fixed=TRUE)
    setNames(sapply(spl, "[", 2), sapply(spl, "[", 1)) })

# check the names of that thing
nam <- unique(unlist(lapply(tab9_list, names)))
stopifnot(length(nam) == length(attrib), all(sort(nam) == sort(attrib)))

# turn into a list with all of the attributes (with NAs as needed)
tab9_df <- lapply(attrib, function(lab) sapply(tab9_list, "[", lab))
names(tab9_df) <- attrib

# add to data as columns
tab <- cbind(tab[,1:8], as.data.frame(tab9_df, stringsAsFactors=FALSE))

# write to database
dbfile <- "mouse_genes.sqlite"
if(file.exists(dbfile)) unlink(dbfile)
db <- dbConnect(SQLite(), dbfile)
dbWriteTable(db, "genes", tab)
dbGetQuery(db, "CREATE INDEX chr_start_stop ON genes (chr, start, stop)")

# add description table
description <- data.frame(description="mouse gene information",
                          source="Mouse Genome Informatics (MGI), Jackson Lab",
                          url=url,
                          date_created=as.character(Sys.Date()),
                          date_source=date_source,
                          genome_build=genome_build,
                          stringsAsFactors=FALSE)
dbWriteTable(db, "description", description, append=TRUE)

dbDisconnect(db)
