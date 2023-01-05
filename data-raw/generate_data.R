glomel <- as.data.frame(readxl::read_excel("data-raw/Glomel.xlsx",sheet = "Feuil1"))
glomel <- glomel[!is.na(glomel$ID), -3]
glomel$Ref <- as.logical(glomel$Ref)
usethis::use_data(glomel, overwrite = TRUE)
