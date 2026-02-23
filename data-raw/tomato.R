## code to prepare `tomato` dataset goes here

## load the tree and data from:
# HIBBINS, Mark S., BREITHAUPT, Lara C. y HAHN, Matthew W.
# Phylogenomic comparative methods: Accurate evolutionary inferences in the presence of gene tree discordance.
# Proceedings of the National Academy of Sciences, 2023, 2022.11.14.516436. (120)

## tree from:
# https://github.com/mhibbins/genetreepruningalg/blob/e2627a711e088844f19c680262e05c433ac81a45/analyses/tomato_analysis/make_tomato_test_input_files.py#L13
tomato_tree <- read.tree(text = "(((((((S.galapagense:60000.0,S.cheesmaniae:60000.0):140000.0,(S.pim1269:100000.0,S.pim1589:100000):100000):1500000.0,(S.neorickii:600000.0,S.arcanum:600000.0):1100000.0):400000.0,(((S.huaylasense:1400000.0,S.peruvianum:1400000.0,S.corneliomulleri:1400000.0):300000.0):1.0,S.chilense:1700000.0):400000.0):300000.0,(S.habrochaites:2000000.0,(S.pen3778:700000.0,S.pen0716:700000):1300000):400000.0):5600000.0):1.0,S.tuberosum:8000000.0);")
# convert from years to coalescent units
tomato_tree[["edge.length"]] <- tomato_tree[["edge.length"]] / 400000
# load the trait
tomato_traits <- read.csv("https://raw.githubusercontent.com/larabreithaupt/seastaR/refs/heads/main/analyses/flower_morphometrics2.csv")

################################################################################
## Format data : keep only one accession number per species

tomato_traits$sp_accession <- paste0(tomato_traits$Sp_ID, "_", tomato_traits$AccessionID)
## accession for each species
aasp <- sapply(unique(tomato_traits$Sp_ID), function(x) summary(factor(subset(tomato_traits, Sp_ID == x)$AccessionID)))
aasp
# when there are several accession per species, we chose the one with the largest number of individual measurements.
# this happens for the following species
aasp[sapply(aasp, function(x) length(x) > 1)]
# tip names vs traits
poptospecies = data.frame(
  population = c( # not included : "S.huaylasense", "S.tuberosum"
    "S.galapagense",  "S.cheesmaniae",  "S.pim1269",         "S.pim1589",
    "S.neorickii",    "S.arcanum",      "S.peruvianum",      "S.corneliomulleri",
    "S.chilense",     "S.habrochaites", "S.pen3778",         "S.pen0716"),
  species_id  = c( # not included : "SCER", "SCHM", "SLYC"
    "SGAL_LA0436",    "SCHE_LA3124",   "SPIM_LA1269",         "SPIM_LA1589",
    "SNEO_LA1321",    "SARC_LA0385",   "SPER_LA0153",         "SCOR_LA0107",
    "SCHL_LA1971",    "SHAB_LA1777",   "SPEN_LA3778",         "SPEN_LA0716")
)

# trim tree to data
tomato_tree_trim <- keep.tip(tomato_tree, poptospecies$population)
tomato_tree_trim$tip.label[match(poptospecies$population, tomato_tree_trim$tip.label)] <- poptospecies$species_id
plot(tomato_tree)
plot(tomato_tree_trim)

# trim data to tree
tomato_traits_trim <- tomato_traits[tomato_traits$sp_accession %in% poptospecies$species_id, ]
tomato_traits_trim <- tomato_traits_trim[do.call(c, sapply(tomato_tree_trim$tip.label, function(x) grep(x, tomato_traits_trim$sp_accession))), ]
tomato_traits_trim$sample_id <- make.unique(tomato_traits_trim$sp_accession)
rownames(tomato_traits_trim) <- tomato_traits_trim$sample_id

tomato_data <- data.frame(Sp_AccessionID = tomato_traits_trim$sp_accession,
                          SampleID = tomato_traits_trim$sample_id,
                          CD = tomato_traits_trim$CD_PlantMean,
                          AL = tomato_traits_trim$AL_PlantMean,
                          SL = tomato_traits_trim$SL_PlantMean
                          )
rownames(tomato_data) <- tomato_data$SampleID

# dataset
tomato <- list(phy = tomato_tree_trim,
               data = tomato_data)

usethis::use_data(tomato, overwrite = TRUE)
