# create a tiny example SQLite database with just the MGI genes that overlap two small regions:
#   chr 2, 2 Mbp interval centered at 97.5 Mbp
#   chr 3, 2 Mbp interval centered at 15.0 Mbp

library(RSQLite)
db <- dbConnect(SQLite(), "mouse_genes.sqlite")
tab <- dbGetQuery(db, paste("SELECT * FROM genes WHERE",
                            "((chr=='2' AND",
                            "((start >= 96500000 AND start <= 98500000) OR",
                            "(stop >=   96500000 AND stop <=  98500000) OR",
                            "(start <=  98500000 AND stop >=  96500000))) OR",
                            "(chr=='3' AND",
                            "((start >= 14000000 AND start <= 16000000) OR",
                            "(stop >=   14000000 AND stop <=  16000000) OR",
                            "(start <=  16000000 AND stop >=  14000000))))",
                            "AND source=='MGI'"))
description <- dbGetQuery(db, "SELECT * FROM description")
dbDisconnect(db)

# write to new database
dbfile <- "../extdata/mouse_genes_small.sqlite"
if(file.exists(dbfile)) unlink(dbfile)
db <- dbConnect(SQLite(), dbfile)
dbWriteTable(db, "genes", tab)
dbGetQuery(db, "CREATE INDEX chr_start_stop ON genes (chr, start, stop)")

# add description table
description[1] <- paste(description[1], "(subset to 2 regions)")
dbWriteTable(db, "description", description, append=TRUE)

dbDisconnect(db)
