## Load data from zip file.
phenoFile <- system.file("extdata", "grainYield_components_BLUEs.zip",
                         package = "statgenGxE")
dropsPheno <- read.csv(unz(description = phenoFile,
                           filename = "2b-GrainYield_components_BLUEs_level.csv"))
## Load genotype meta data.
genoMeta <- read.csv(system.file("extdata", "8-Info_Maize_variety.csv",
                                 package = "statgenGxE"))
## Rename genetic_group geneticGroup for consistency.
colnames(genoMeta)[colnames(genoMeta) == "genetic_group"] <- "geneticGroup"

## Restrict to 10 relevant environments.
exps <- c("Cam12R", "Cra12R", "Gai12W", "Kar12W", "Kar13R", "Kar13W",
          "Mar13R", "Mur13R", "Mur13W", "Ner12R")
dropsPheno <- dropsPheno[dropsPheno[["Experiment"]] %in% exps, ]

## Add year column.
dropsPheno[["year"]] <- paste0("20", substring(dropsPheno[["Experiment"]],
                                               first = 4, last = 5))
## Add location column.
dropsPheno[["loc"]] <- substring(dropsPheno[["Experiment"]],
                                 first = 1, last = 3)

## Add scenario columns.
scenario <- data.frame(Experiment = exps,
                       scenarioWater = c("WD", "WD", "WW", "WW", "WD",
                                         "WW", "WD", "WW", "WW", "WD"),
                       scenarioTemp = c("Hot", "Hot", "Cool", "Cool",
                                        "Hot(Day)", "Hot(Day)", "Hot(Day)",
                                        "Hot", "Hot", "Hot(Day)"))
dropsPheno <- merge(dropsPheno, scenario)

dropsPheno[["scenarioFull"]] <- interaction(dropsPheno[c("scenarioWater",
                                                          "scenarioTemp")],
                                            drop = TRUE)

## Add genetic groups.
dropsPheno <- merge(dropsPheno, genoMeta[c("Variety_ID", "geneticGroup")])

dropsPheno <- droplevels(dropsPheno)

usethis::use_data(dropsPheno, overwrite = TRUE)
