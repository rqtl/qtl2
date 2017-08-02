# create sqlite database with mouse genes, from MGI
#
# source (dated 2016-02-24)
#      ftp://ftp.jax.org/SNPtools/genes/MGI.sorted.txt.gz
#      ftp://ftp.jax.org/SNPtools/genes/MGI.sorted.txt.gz.tbi
#  (.tbi file not used here)

library(RSQLite)

### download files
site <- "ftp://ftp.jax.org/SNPtools/genes"
file <- "MGI.sorted.txt.gz"
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
fields <- c("seqid", "source", "type", "start", "stop", "score", "strand",
            "phase", "ID", "Name", "Parent", "Dbxref", "mgiName", "bioType")

# load file
tab <- data.table::fread(file, data.table=FALSE,
                         colClasses=rep(c("character", "numeric", "character"), c(3,2,9)))
colnames(tab) <- fields

if(remove_tmpfile) { # remove temporary file
    unlink(tmpfile)
}

# remove things not on chr 1-19, "X", "Y", "MT"
tab <- tab[tab$seqid %in% c(1:19, "X", "Y", "MT"),]

# score -> numeric
tab$score[tab$score=="."] <- NA
tab$score <- as.numeric(tab$score)
tab$phase[tab$phase=="."] <- NA
tab$phase <- as.numeric(tab$phase)

# strand: . -> missing
tab$strand[tab$strand=="."] <- NA

# seqid -> chr
names(tab)[1] <- "chr"

# write to database
dbfile <- "mouse_genes.sqlite"
if(file.exists(dbfile)) unlink(dbfile)
db <- dbConnect(SQLite(), dbfile)
dbWriteTable(db, "genes", tab)
dbGetQuery(db, "CREATE INDEX chr_start ON genes (chr, start)")
dbGetQuery(db, "CREATE INDEX chr_stop ON genes (chr, stop)")
dbGetQuery(db, "CREATE INDEX chr_start_stop ON genes (chr, start, stop)")

# add description table
description <- data.frame(description="mouse gene information",
                          source="Mouse Genome Informatics (MGI), Jackson Lab",
                          url=url,
                          date_created=as.character(Sys.Date()),
                          date_source="2016-02-24",
                          genome_build="GRCm38/mm10")
dbWriteTable(db, "description", description, append=TRUE)

dbDisconnect(db)
