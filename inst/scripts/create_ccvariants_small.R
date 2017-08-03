# create a tiny example SQLite database with just the variants that overlap two small regions:
#   chr 2, 2 Mbp interval centered at 97.5 Mbp
#   chr 3, 2 Mbp interval centered at 15.0 Mbp

library(RSQLite)
db <- dbConnect(SQLite(), "cc_variants.sqlite")

# description
description <- dbGetQuery(db, "SELECT * FROM description")
description$description <- paste(description$description, "(subset to 2 regions)")

tab <- dbGetQuery(db, paste("SELECT * FROM variants WHERE",
                            "(chr=='2' AND pos >= 96500000 AND pos <= 98500000) OR",
                            "(chr=='3' AND pos >= 14000000 AND pos <= 16000000)"))
dbDisconnect(db)

# write to new database
dbfile <- "../extdata/cc_variants_small.sqlite"
if(file.exists(dbfile)) unlink(dbfile)
db <- dbConnect(SQLite(), dbfile)
dbWriteTable(db, "variants", tab)
dbGetQuery(db, "CREATE INDEX chr_pos ON variants (chr, pos)")

# add description table
dbWriteTable(db, "description", description, append=TRUE)

dbDisconnect(db)
