
## For writing plink's file formats

fam <- function(df, trait) data.frame(fid = df$animal_id,
                                      iid = df$animal_id,
                                      father = 0,
                                      mother = 0,
                                      sex = 2,
                                      as.data.frame(df[, trait]))

write_plink <- function(x, filename) write.table(x,
                                                 file = filename,
                                                 quote = FALSE,
                                                 row.names = FALSE,
                                                 col.names = FALSE)