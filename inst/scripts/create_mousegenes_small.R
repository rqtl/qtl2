# create a tiny example SQLite database with just the MGI genes that overlap two small regions:
#   chr 2, 2 Mbp interval centered at 97.5 Mbp
#   chr 3, 2 Mbp interval centered at 15.0 Mbp

library(RSQLite)
db <- dbConnect(SQLite(), "mouse_genes.sqlite")
tab <- dbGetQuery(db, paste("SELECT * FROM genes WHERE",
                            "((chr=='2' AND",
                            "((start_Mbp >= 97.5-1 AND start_Mbp <= 97.5-1) OR",
                            "(stop_Mbp >= 97.5-1 AND stop_Mbp <= 97.5-1) OR",
                            "(start_Mbp <= 97.5+1 AND stop_Mbp >= 97.5-1))) OR",
                            "(chr=='3' AND",
                            "((start_Mbp >= 15.0-1 AND start_Mbp <= 15.0-1) OR",
                            "(stop_Mbp >= 15.0-1 AND stop_Mbp <= 15.0-1) OR",
                            "(start_Mbp <= 15.0+1 AND stop_Mbp >= 15.0-1))))",
                            "AND source=='MGI'"))
dbDisconnect(db)

# write to new database
dbfile <- "../extdata/mouse_genes_small.sqlite"
if(file.exists(dbfile)) unlink(dbfile)
db <- dbConnect(SQLite(), dbfile)
dbWriteTable(db, "genes", tab)
dbGetQuery(db, "CREATE INDEX chr_start_Mbp ON genes (chr, start_Mbp)")
dbGetQuery(db, "CREATE INDEX chr_stop_Mbp ON genes (chr, stop_Mbp)")
dbGetQuery(db, "CREATE INDEX chr_start_stop_Mbp ON genes (chr, start_Mbp, stop_Mbp)")

# add description table
description <- data.frame(description="mouse gene information (subset to 2 regions)",
                          source="Mouse Genome Informatics (MGI), Jackson Lab",
                          url=url,
                          date_created=as.character(Sys.Date()),
                          date_source="2016-02-24",
                          genome_build="GRCm38/mm10")
dbWriteTable(db, "description", description, append=TRUE)

dbDisconnect(db)
