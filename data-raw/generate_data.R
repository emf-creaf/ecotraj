glomel <- as.data.frame(readxl::read_excel("data-raw/Glomel.xlsx",sheet = "Feuil1"))
glomel <- glomel[!is.na(glomel$ID), -c(3,4)]
glomel$Ref <- as.logical(glomel$Ref)
usethis::use_data(glomel, overwrite = TRUE)

load("data-raw/glenan")
usethis::use_data(glenan, overwrite = TRUE)
