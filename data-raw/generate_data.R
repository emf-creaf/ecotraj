glomel <- as.data.frame(readxl::read_excel("data-raw/Glomel.xlsx",sheet = "Feuil1"))
glomel <- glomel[!is.na(glomel$ID), -3]
glomel$Ref <- as.logical(glomel$Ref)
usethis::use_data(glomel, overwrite = TRUE)

# glenan <- as.data.frame(readxl::read_excel("data-raw/Glenan.xlsx"))
# glenan <- glenan[seq(2,104, by=2), ]
# row.names(glenan) <-NULL # Reset row names
# usethis::use_data(glenan, overwrite = TRUE)
