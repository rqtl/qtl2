# create sqlite database with mouse genes, from MGI
#
# source (dated 2017-08-03):
#      http://www.informatics.jax.org/downloads/mgigff/
#          MGI.20170803.gff3.gz
#          MGI_GFF_Spec.docx (annotations, e.g. column names)

library(RSQLite)

### download files
site <- "http://www.informatics.jax.org/downloads/mgigff"
file <- "MGI.gff3.gz"
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
attrib <- c("ID", "Name", "Parent", "Dbxref", "mgiName", "bioType")

# read data
tab <- read.delim(file, sep="\t", header=FALSE, comment.char="#",
                  na.strings=".", stringsAsFactors=FALSE,
                  colClasses=rep(rep(c("character", "numeric"), 3), c(3,3,1,1,1,0)))
colnames(tab) <- fields

# remove things not on chr 1-19, "X", "Y", "MT"
tab <- tab[tab$chr %in% c(1:19, "X", "Y", "MT"),]

# split 9th column at ';' then split at '=' and use first bit as key and second bit as value
tab9_spl <- strsplit(tab[,9], ";")
tab9_list <- lapply(tab9_spl, function(a) {
    spl <- strsplit(a, "=", fixed=TRUE)
    setNames(sapply(spl, "[", 2), sapply(spl, "[", 1)) })

# turn into a list with all of the attributes (with NAs as needed)
tab9_df <- lapply(attrib, function(lab) sapply(tab9_list, "[", lab))
names(tab9_df) <- attrib

# add to data as columns
tab <- cbind(tab[,1:8], as.data.frame(tab9_df, stringsAsFactors=FALSE))

# remove temporary file
if(remove_tmpfile) {
    unlink(tmpfile)
}

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
                          date_source="2017-08-03",
                          genome_build="GRCm38/mm10")
dbWriteTable(db, "description", description, append=TRUE)

dbDisconnect(db)
