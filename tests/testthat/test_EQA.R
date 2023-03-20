library(ecotraj)
data(glomel)

# Extract compositional data matrix
glomel_comp <- as.matrix(glomel[,!(names(glomel) %in% c("ID", "Ref", "Complementary"))])
rownames(glomel_comp) <- glomel$ID

# Calculate Bray-curtis distance matrix 
glomel_bc <- vegan::vegdist(glomel_comp, method = "bray")

# Define reference envelope by observation ID
glomel_env <- glomel$ID[glomel$Ref]

test_that("EQA can be performed",{
    # Assess quality with respect to reference envelope
  expect_s3_class(compareToStateEnvelope(glomel_bc, glomel_env), "data.frame")
})