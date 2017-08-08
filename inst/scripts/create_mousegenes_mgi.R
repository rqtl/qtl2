# create a SQLite database just the MGI genes

library(RSQLite)
db <- dbConnect(SQLite(), "mouse_genes.sqlite")
tab <- dbGetQuery(db, "SELECT * FROM genes WHERE source='MGI'")
dbDisconnect(db)

# write to new database
dbfile <- "mouse_genes_mgi.sqlite"
if(file.exists(dbfile)) unlink(dbfile)
db <- dbConnect(SQLite(), dbfile)
dbWriteTable(db, "genes", tab)
dbGetQuery(db, "CREATE INDEX chr_start_stop ON genes (chr, start, stop)")

# add description table
description <- data.frame(description="mouse gene information (subset to source='MGI' records)",
                          source="Mouse Genome Informatics (MGI), Jackson Lab",
                          url=url,
                          date_created=as.character(Sys.Date()),
                          date_source="2016-02-24",
                          genome_build="GRCm38/mm10")
dbWriteTable(db, "description", description, append=TRUE)

dbDisconnect(db)
