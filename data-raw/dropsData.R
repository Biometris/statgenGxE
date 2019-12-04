## Load data from zip file.
phenoFile <- system.file("extdata", "grainYield_components_BLUEs.zip",
                         package = "statgenGWAS")
dropsPheno <- read.csv(unz(description = phenoFile,
                           filename = "2b-GrainYield_components_BLUEs_level.csv"))
## Add year column to demonstrate plot options.
dropsPheno[["year"]] <- paste0("20", substring(dropsPheno[["Experiment"]],
                                               first = 4, last = 5))
## Remove observations from 2011 from dropsPheno.
dropsPheno <- dropsPheno[dropsPheno[["year"]] != "2011", ]
dropsPheno <- droplevels(dropsPheno)

usethis::use_data(dropsPheno, overwrite = TRUE)
