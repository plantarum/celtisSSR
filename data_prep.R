library(polysat)
library(dplyr)
library(adegenet)
library(poppr)
library(tidyr)
library(ggplot2)
library(cowplot)
library(gridGraphics)
library(hierfstat)
library(raster)
library(magrittr)
library(sf)

##############
## SSR Data ##
##############

message("Importing ssr data")
## Import and correct raw data:
cssr <- read.GeneMapper("data/cssr2018-09-06.txt") #reads the file

#################
## Corrections ##
#################

## missed the following peak in Geneious:
Genotype(cssr, sample = "ARCC3", locus = "loc31") <- 453

## lowest allele is smaller, looks like a stutter - shouldn't have 4 alleles for a triploid!!
Genotype(cssr, sample = "SCCC5", locus = "loc28") <- c(200, 206, 212)

## largest allele is smaller than the first three, remove it
Genotype(cssr, sample = "MOSJ11", locus = "loc20") <- c(241, 245, 249)

## we scored a peak at 407, it's a bit small, and this is a diploid 
Genotype(cssr, sample = "WB3", locus = "loc24") <- c(395, 437)

## we scored a peak at 407, it's a bit small, and this is a diploid
Genotype(cssr, sample = "OHCQ2", locus = "loc24") <- c(395, 437) 

## took out third peak, too small and should be diploid
Genotype(cssr, sample = "OHES3", locus = "loc24") <- c(395, 407)

## took out third peak, too small and should be diploid
Genotype(cssr, sample = "ILWB22", locus = "loc38") <- c(346, 352)

## took out two peaks, too small and should be diploid
Genotype(cssr, sample = "INEM12", locus = "loc38") <- c(328)

## took out third peak, two small and should be diploid
Genotype(cssr, sample = "ALCC10", locus = "loc38") <- c(346, 358)

Usatnts(cssr) <- c(3,3,6,6,6,6,6,6,6,6) #assigns the repeat number of nucleotides

cssr <- deleteLoci(cssr, c("loc39")) # this primer pair was messed up in
                                     # the lab, untrustworthy

## There were a lot of diploids with three alleles for locus 24, so
## something fishy is going on with that one. We'll drop it for now.
cssr <- deleteLoci(cssr, c("loc24"))

## Clean out bad samples - anyone with more than 3 primer pairs missing
misses <- find.missing.gen(cssr)
misses <- misses %>% group_by(Sample) %>% summarise(total = n()) %>% as.data.frame()
misses <- misses[order(misses$total),]

cssr <- deleteSamples(cssr, misses[misses$total > 2, "Sample"])

# deleted samples that had four peaks in SSR data even though they are
# triploid according to flow (at loc16,20) 
cssr <- deleteSamples(cssr, samples = c("KYHR5", "ILWB5"))

###################################
## Import metadata into popTable ##
###################################

###############################################################################
## NB: popTable contains species determinations. These are unreliable!! They ##
## are based on our field ID, which, especially for small plants, were       ##
## guesses. We will assign species names based on ploidy and SSR clusters    ##
## later on in this script.                                                  ##
###############################################################################

message("Importing metadata into popTable")
popTable <- read.csv("data/popTable.csv", as.is = FALSE)
popTable$sample <- as.character(popTable$sample)
row.names(popTable) <- as.character(popTable$sample)
coordinates(popTable) <- ~Longitude+Latitude
crs(popTable) <- CRS("+proj=longlat +datum=WGS84")

## These two samples were categorized as triploid. However, the flow data
## was very low quality, and on review doesn't look trustworthy
popTable["VARC2", "ploidy"] <- "unknown"
popTable["VARC3", "ploidy"] <- "unknown"

## Correct the locations for Pinery Provincial Park
popTable[popTable$population == "PPP", "Latitude"] <- 43.25767
popTable[popTable$population == "PPP", "Longitude"] <- -81.835987

## Sideling Hill location correction:
popTable[popTable$population == "PASH", "Latitude"] <- 39.732734
popTable[popTable$population == "PASH", "Longitude"] <- -78.347700

## TWS was collected at Blue Licks Battlefield State Park:
popTable[popTable$population == "Tws", "Latitude"] <- 38.451706
popTable[popTable$population == "Tws", "Longitude"] <- -83.993636

popTable["ILHE1", "population"] <- "ILHE"

popTable$population <- factor(popTable$population)

## Add a column to popTable that has 2 for diploid, 3 for triploids:
popTable$ploidy <- as.character(popTable$ploidy)
popTable$ploidy[popTable$ploidy == "diploid"] <- 2
popTable$ploidy[popTable$ploidy == "triploid"] <- 3
popTable$ploidy[popTable$ploidy == "tetraploid"] <- 4
popTable$ploidy[popTable$ploidy == "unknown"] <- NA

popTable$ploidy <- as.numeric(popTable$ploidy)

## combine the three different Pt Pelee 'populations'
levels(popTable$population) <- c(levels(popTable$population), "ONPE")
popTable$population[popTable$population %in% c("COHA", "VCN", "WB")] <- "ONPE"
popTable$population <- factor(popTable$population)

################################################################################
## Various different visualizations used different regions, areas, species,   ##
## populations, etc. I dealt with this by adding columns to popTable. Not all ##
## of this turned out to be interesting, but I leave it here in case we want  ##
## to follow up on any of this.                                               ##
################################################################################

## Assign regions and areas to populations:
popTable$region  <-  "MW"

popTable$region[popTable$population %in%
                c("ILGG", "ILHE", "ILLG", "ILWB")] <- "IL" 
popTable$region[popTable$population %in%
                c("ONPE", "GC", "PtA", "Landsdale", "PPP")] <- "ON" 
popTable$region[popTable$population %in%
                c("ALAC", "ALCC", "ALHR", "SCCC", "SCSF")] <- "SE" 
## popTable$region[popTable$population %in%
##                 c("KYHR", "Tws", "VAWR")] <- "KY"
popTable$region[popTable$population %in%
                c("PASH", "PABE")] <- "PA"
popTable$region[popTable$population %in%
                c("OHCQ", "OHES", "MIWA")] <- "GL"
popTable$region[(popTable$population == "TXFW") &
                (popTable$ploidy == 3)]  <- "C. reticulata"
popTable$region[popTable$population %in%
                c("MSCS", "MSSM", "MSWD")] <- "MS"
popTable$region[popTable$population %in%
                c("VADS", "VARC")] <- "VA"
## popTable$region[popTable$population %in%
##                 c("SCCC", "SCSF")] <- "SC"

popTable$area <- "S"
popTable$area[popTable$region %in% c("ON", "GL", "PA")] <- "N"
popTable$area[popTable$region %in% c("C. reticulata")] <- "C. reticulata"

##################################
## Import STRUCTURE assignments ##
##################################
message("Importing STRUCTURE assignments")
##############################################################
## cluster/admixture assignments as generated by STRUCTURE: ##
##############################################################
dipAssign <- read.table("data/dipAssign.csv")
dipAssign$population  <- as.factor(dipAssign$population) # because why not?
dipAssign$sample <- as.character(dipAssign$sample)
dipAssign$location <- droplevels(popTable[dipAssign$sample, ]$population)


popOrder <- dipAssign %>% mutate(one = pop1) %>%
  gather(key = pop, value = prop, pop1, pop2) %>%
  group_by(location) %>% summarize(mean = mean(one)) %>%
  arrange(mean)

dipAssign$location <- factor(dipAssign$location, levels = popOrder$location)

popTable$species <- NA ## original values were guesses
popTable[dipAssign$sample, "species"] <- as.character(dipAssign$population)
popTable[dipAssign$sample, "area"] <- as.character(dipAssign$population)
popTable[dipAssign$sample, "region"] <- as.character(dipAssign$population)
popTable[which(popTable$ploidy == 3), "species"] <- "pumila"
popTable[which(popTable$ploidy == 4), "species"] <- "tetraploid"

############################
## Apply metadata to cssr ##
############################

message("applying metadata to cssr")
Description(cssr) <- "Dataset for Celtis SSRs"   #creates a title for data set
PopNames(cssr) <- levels(popTable$population)    #designates names for the populations
PopInfo(cssr) <-  as.numeric(popTable[Samples(cssr), ]$population)

cssr <- reformatPloidies(cssr, output = "sample")

## subset the popTable to only include samples that had usable SSR
## fingerprints: 
popTableInCSSR <- popTable[popTable$sample %in% Samples(cssr), ]
Ploidies(cssr)[as.character(popTableInCSSR$sample)] <- popTableInCSSR$ploidy

## pull out diploid samples:
diploids <- Samples(cssr, ploidies = 2)
diploids <- cssr[diploids]
PopNames(diploids) <- levels(dipAssign$population)

dipAssignInCSSR <- dipAssign[dipAssign$sample %in% Samples(cssr), ]
PopInfo(diploids)[dipAssignInCSSR$sample]  <- as.numeric(dipAssignInCSSR$population)

## making a subset of samples from just the triploids in the dataset
triploids <- Samples(cssr, ploidies = 3)
triploids <- cssr[triploids]

PopNames(triploids) <- levels(factor(c("South", "North", "C. reticulata")))
PopInfo(triploids)  <- factor(popTable[Samples(triploids), ]$area)

allPloidy <- merge(diploids, triploids) 

tetraploids <- Samples(cssr, ploidies = 4)
tetraploids <- cssr[tetraploids]
PopInfo(tetraploids) <- 1
PopNames(tetraploids) <- "tetraploid"

## only tetraploids
justTetra <- tetraploids

## all data, including tetraploids:
tetraploids <- merge(allPloidy, tetraploids)

########################
## Export to adegenet ##
########################
message("Exporting to adegenet")

triploidAG <- gendata.to.genind(triploids)
triploidAGStrata <-
  data.frame(species = pop(triploidAG),
             population = popTable[rownames(triploidAG@tab), "population"])
strata(triploidAG) <- triploidAGStrata

#########################
## Export to STRUCTURE ##
#########################

#################################################################################
## The following mash up exports the ~diploids~ genambig object to a STRUCTURE ##
## file, which then allows us to import it into a genind object, ~diploidsGI~. ##
## That seems rather round about, but I couldn't find a direct genambig-genind ##
## conversion.                                                                 ##
#################################################################################

message("Exporting to STRUCTURE")
write.Structure(diploids, ploidy = 2, file = "data/diploidStructure.stru")
tmpFile <- readLines("data/diploidStructure.stru")
tmpFile[1] <- gsub("rowlabel\t", "", tmpFile[1])
tmpFile[1] <- gsub("PopInfo\t", "", tmpFile[1])
tmpFile[2] <- gsub("missing\t", "", tmpFile[2])
writeLines(text = tmpFile, con = "data/diploidStructure.stru")

##############################################
## Import from STRUCTURE to adegenet format ##
##############################################
message("Importing from STRUCTURE")
diploidGI <- read.structure("data/diploidStructure.stru",
                            n.ind = length(Samples(diploids)),
                            n.loc = length(Loci(diploids)),
                            onerowperind = FALSE, col.lab = 1, col.pop = 2,
                            ask = FALSE)
popNames(diploidGI) <- c("occidentalis", "laevigata", "mix")
locNames(diploidGI) <- Loci(diploids)
diploidStrata <- data.frame(species = pop(diploidGI),
                            population = popTable[rownames(diploidGI@tab), ]$population)
strata(diploidGI) <- diploidStrata

diploidPop <- genind2genpop(diploidGI)

diploidOL <- diploidGI[strata(diploidGI)$species != "mix"]
setPop(diploidOL) <- ~species

diploidSep <- seppop(diploidOL)

setPop(diploidSep[[1]]) <- ~population
setPop(diploidSep[[2]]) <- ~population

pop(diploidSep[[1]]) <- popTable[rownames(diploidSep[[1]]@tab), ]$population
pop(diploidSep[[2]]) <- popTable[rownames(diploidSep[[2]]@tab), ]$population

diploidHybrids <- diploidGI[strata(diploidGI)$species == "mix"]
setPop(diploidHybrids) <- ~population
## diploidSep is a list containing two GI objects, one for each diploid
## species. They are divided into populations.

## Diversity Summaries

bsLaev <- basic.stats(diploidSep$laevigata)
LaevStats <- data.frame(species = "laevigata",
                       N = NA,
                       MLG = NA,
                       Ho = round(colMeans(bsLaev$Ho, na.rm = TRUE), 3),
                       Hs = round(colMeans(bsLaev$Hs, na.rm = TRUE), 3), 
                       Fis = round(colMeans(bsLaev$Fis, na.rm = TRUE), 3))
nLaev <- data.frame(n = table(pop(diploidSep$laevigata)))
nLaev$n.Var1 <- as.character(nLaev$n.Var1)
LaevStats[nLaev[, 1], "N"] <- nLaev[, 2]
LaevStats$MLG  <- data.frame(MLG =
                              apply(mlg.table(diploidSep$laevigata, plot = FALSE),
                                    MAR = 1, function(x) sum(as.logical(x))))

bsOcci <- basic.stats(diploidSep$occidentalis)
OcciStats <- data.frame(species = "occidentalis",
                       N = NA,
                       MLG = NA,
                       Ho =  round(colMeans(bsOcci$Ho, na.rm = TRUE), 3),
                       Hs =  round(colMeans(bsOcci$Hs, na.rm = TRUE), 3),
                       Fis = round(colMeans(bsOcci$Fis, na.rm = TRUE), 3))
nOcci <- data.frame(n = table(pop(diploidSep$occidentalis)))
nOcci$n.Var1 <- as.character(nOcci$n.Var1)
OcciStats[nOcci[, 1], "N"] <- nOcci[, 2]
OcciStats$MLG  <- data.frame(MLG =
                              apply(mlg.table(diploidSep$occidentalis,
                                              plot = FALSE),
                                    MAR = 1, function(x) sum(as.logical(x))))

bsHybrids <- basic.stats(diploidHybrids)
HybridStats <- data.frame(species = "hybrid",
                       N = NA,
                       MLG = NA,
                       Ho =  round(colMeans(bsHybrids$Ho, na.rm = TRUE), 3),
                       Hs =  round(colMeans(bsHybrids$Hs, na.rm = TRUE), 3),
                       Fis = round(colMeans(bsHybrids$Fis, na.rm = TRUE), 3))
nHybrids <- data.frame(n = table(pop(diploidHybrids)))
nHybrids$n.Var1 <- as.character(nHybrids$n.Var1)
HybridStats[nHybrids[, 1], "N"] <- nHybrids[, 2]
HybridStats$MLG  <- data.frame(MLG =
                              apply(mlg.table(diploidHybrids, plot = FALSE),
                                    MAR = 1, function(x) sum(as.logical(x))))

## basic.stats:
##  Hs=\tilde{n}/(\tilde{n}-1)[1-sum_i\bar{p_i^2}-Ho/2\tilde{n}],

##      where \tilde{n}=np/sum_k 1/n_k and \bar{p_i^2}=sum_k p_{ki}^2/np

## adegenet:

##   Let _m(k)_ be the number of alleles of locus _k_, with a total of
##      _K_ loci. We note f_i the allele frequency of allele _i_ in a
##      given population. Then, Hs is given for a given population by:
##      \frac{1}{K} sum_{k=1}^K (1 - sum_{i=1}^{m(k)} f_i^2)


## Allele frequencies
allelefreq <- simpleFreq(allPloidy)
occiAlleles <- allelefreq["occidentalis", -1] > 0 
laevAlleles <- allelefreq["laevigata", -1] > 0 
northAlleles <- allelefreq["North", -1] > 0  
southAlleles <- allelefreq["South", -1] > 0  
reticAlleles <- allelefreq["C. reticulata", -1] > 0  

occiPA <- occiAlleles & (! laevAlleles)
laevPA <- laevAlleles & (! occiAlleles)

dipStats <- data.frame(Ap = numeric(), Ho = numeric(), Hs = numeric(),
                      Fis = numeric(), Fst = numeric())
dipStats["occidentalis", ] <- c(sum(occiPA),
                               bsOcci$overall[c("Ho", "Hs", "Fis", "Fst")])
dipStats["laevigata", ] <- c(sum(laevPA),
                            bsLaev$overall[c("Ho", "Hs", "Fis", "Fst")])

northOcci <- sum(northAlleles & occiPA)
northLaev <- sum(northAlleles & laevPA)

southOcci <- sum(southAlleles & occiPA)
southLaev <- sum(southAlleles & laevPA)

pwFst <- calcPopDiff(allelefreq[row.names(allelefreq) != "mix", ],
                    metric = "Fst")

pwFst <- pwFst[c("laevigata", "occidentalis", "North", "South",
                "C. reticulata"),
              c("laevigata", "occidentalis", "North", "South",
                "C. reticulata")] 

simGst <- calcPopDiff(allelefreq, metric = "Gst")

########################
## MLG/Clone Analysis ##
########################

tripDist <- meandistance.matrix(triploids,
                               distmetric = Bruvo2.distance,
                               add = TRUE, loss = FALSE)
cloneDF <- data.frame(sample = Samples(triploids),
                      population =
                        PopNames(triploids)[PopInfo(triploids)], 
                      mlg = assignClones(tripDist,
                                         threshold = 0.1))

cloneTable <- table(cloneDF[, c("population", "mlg")])

clonePA <- cloneTable > 0

popClones <- rowSums(clonePA)

clonePops <- colSums(clonePA)

PumilaMLG <- data.frame(species = "pumila",
                       N = rowSums(cloneTable),
                       MLG = as.numeric(popClones))


###########################
## Triploids to adegenet ##
###########################
## This doesn't work, as read.structure requires haploid or diploid
## samples.
## 
## triploidGI <-
##   read.structure("triploidStructure.stru",
##                  n.ind = length(Samples(triploids)),  
##                  n.loc = length(Loci(triploids)),
##                  onerowperind = FALSE, col.lab = 1, col.pop = 2,
##                  ask = FALSE)
## popNames(triploidGI) <- c("tenuifolia")
## locNames(triploidGI) <- Loci(triploids)
## triploidStrata <-
##   data.frame(species = pop(triploidGI),
##              population = popTable[rownames(triploidGI@tab),
##                                    "population"]) 
## strata(triploidGI) <- triploidStrata

################
## Morphology ##
################

###############################################
## This is not used in the manuscript (yet?) ##
###############################################

## message("Morphology")
## morph <- read.csv("data/Field Measurements.csv", row.names = 1)
## ## cleanup leading spaces from row-names:
## row.names(morph) <- gsub("^ +", "", row.names(morph))
## morph$species <- NA

## popTableInMorph <- popTable[popTable$sample %in% row.names(morph), ]
## morph[popTableInMorph$sample, "species"] <- popTableInMorph$species

## morph$densRL <- morph$roundlower / (morph$length/2)
## morph$densSL <- morph$sharplower / (morph$length/2)
## morph$densRU <- morph$roundupper / (morph$length/2)
## morph$densSU <- morph$sharpupper / (morph$length/2)

## morphComp <- morph[complete.cases(morph), ]
## morphDip <- morphComp[morphComp$species %in%
##                       c("laevigata", "occidentalis", "mix"), ]
## morphDip$species <- factor(morphDip$species)

## morphOL <- morphDip[morphDip$species != "mix", ]
## morphOL$species <- factor(morphOL$species)

#############
## Mapping ##
#############
message("prepping maps")

us <- getData("GADM", country = "USA", level = 1, path = "./data/maps/",
             type = "sf")
canada <- getData("GADM", country = "CAN", level = 1, path = "./data/maps",
                 type = "sf")
mex <- getData("GADM", country = "MEX", level = 1, path = "./data/maps",
              type = "sf")

na <- rbind(us, canada, mex)
suppressWarnings(na.simp <- st_simplify(na, dTolerance = 0.01))
laea = CRS("+proj=laea +lat_0=30 +lon_0=-95")
na.la <- st_transform(na.simp, laea)

greatlakes <- st_read("data/maps/greatlakes.shp")
gl.proj <- st_transform(greatlakes, laea)

popTableProj <- spTransform(popTable, laea)

flow <- read.csv("data/flow.csv")
flow$pg <- gsub("#VALUE!", NA, flow$pg)
flow$pg <- as.numeric(flow$pg)
flowSum <- list(trip = list(), dip = list(), tet = list())
flowSum$trip <- list(mean = mean(subset(flow, ploidy == "triploid")$pg),
                    sd = sd(subset(flow, ploidy == "triploid")$pg))
flowSum$dip <- list(mean = mean(subset(flow, ploidy == "diploid")$pg),
                    sd = sd(subset(flow, ploidy == "diploid")$pg))
flowSum$tet <- list(mean = mean(subset(flow, ploidy == "tetraploid")$pg),
                    sd = sd(subset(flow, ploidy == "tetraploid")$pg))

message("finished processing data")
message("dipStats has ", nrow(dipStats), " rows")
