## Load data from zip file.
phenoFile <- system.file("extdata", "grainYield_components_BLUEs.zip",
                         package = "statgenGWAS")
dropsPheno <- read.csv(unz(description = phenoFile,
                           filename = "2b-GrainYield_components_BLUEs_level.csv"))
## Remove observations from 2011 from dropsPheno.
dropsPheno <- dropsPheno[substring(dropsPheno[["Experiment"]],
                                   first = 4, last = 5) != "11" &
                           dropsPheno[["Experiment"]] != "Gra13W", ]
dropsPheno <- droplevels(dropsPheno)

usethis::use_data(dropsPheno, overwrite = TRUE)
