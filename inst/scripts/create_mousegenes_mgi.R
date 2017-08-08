# create a SQLite database just the MGI genes

library(RSQLite)
db <- dbConnect(SQLite(), "mouse_genes.sqlite")
tab <- dbGetQuery(db, "SELECT * FROM genes WHERE source='MGI'")
description <- dbGetQuery(db, "SELECT * FROM description")
dbDisconnect(db)

# write to new database
dbfile <- "mouse_genes_mgi.sqlite"
if(file.exists(dbfile)) unlink(dbfile)
db <- dbConnect(SQLite(), dbfile)
dbWriteTable(db, "genes", tab)
dbGetQuery(db, "CREATE INDEX chr_start_stop ON genes (chr, start, stop)")

# add description table
description[1] <- paste(description[1], "(subset to source='MGI' records)")
dbWriteTable(db, "description", description, append=TRUE)

dbDisconnect(db)
