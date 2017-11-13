# create a SQLite database just the MGI genes

library(RSQLite)
db <- dbConnect(SQLite(), "mouse_genes.sqlite")
tab <- dbGetQuery(db, "SELECT * FROM genes WHERE source='MGI'")
description <- dbGetQuery(db, "SELECT * FROM description")
dbDisconnect(db)

# remove duplicates (they have multiple Dbxref but are otherwise identical)
#    ... pick a random one of each
tab_name <- table(tab$Name)
dup <- names(tab_name)[tab_name > 1]
omit <- NULL
for(i in dup) {
    z <- which(tab$Name==i)
    dbxref <- unique(tab$Dbxref[z])
    keep <- sample(length(z), 1)
    tab$Dbxref[z[keep]] <- paste(dbxref, collapse=",")
    omit <- c(omit, z[-keep])
}
tab <- tab[-omit,]

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
