## library(polysat)
## library(dplyr)
## library(adegenet)
## library(poppr)
## library(tidyr)
## library(ggplot2)
## library(cowplot)
## library(gridGraphics)
## library(hierfstat)
## library(raster)
## library(magrittr)
## library(sf)

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

#########################
## Full Structure Data ##
#########################

message("all STRUCTURE assignments")
## 2 clusters
allAssign2 <- read.csv("structure-all/all2_1.csv", header = FALSE)
allAssign2 <- allAssign2[, -1]
colnames(allAssign2) <- c("sample", "missing", "population",
                         "pop1", "pop2")
rownames(allAssign2) <- allAssign2$sample

allAssign2$ploidy <-
  as.data.frame(popTable)[rownames(allAssign2), "ploidy"]
allAssign2$species <-
  as.data.frame(popTable)[rownames(allAssign2), "species"]
allAssign2$population <-
  as.data.frame(popTable)[rownames(allAssign2), "population"]
allAssign2$popPloid <- paste(allAssign2$population,
                            allAssign2$ploidy, sep = "")
allAssign2$region <-
  as.data.frame(popTable)[rownames(allAssign2), "area"]

allAssign2[c("SCCC5", "MOSJ11"), "region" ] <- "tetraploid"

allAssign2$sortA <- paste(allAssign2$ploidy, allAssign2$species,
                         allAssign2$region,
                         allAssign2$population, sep = "") 

allAssign2 <- allAssign2[order(allAssign2$sortA), ]

xlabels <- aggregate(1:nrow(allAssign2),
                    by = list(allAssign2[, "region"]),
                    FUN = mean)
xlabels[, 1] <- c(expression(italic("Cr")),
                 expression(italic("C. laevigata")),
                 expression("hybrid"),
                 expression("Northern Triploids"),
                 expression(italic("C. occidentalis")),
                 expression("Southern Triploids"),
                 expression("4"))

sampleEdges <- aggregate(1:nrow(allAssign2),
                        by = list(allAssign2[, "region"]), 
                        FUN = max)

## 3 Clusters
allAssign3 <- read.csv("structure-all/all3_1.csv", header = FALSE)
allAssign3 <- allAssign3[, -1]
colnames(allAssign3) <- c("sample", "missing", "population",
                         "pop1", "pop2", "pop3")
rownames(allAssign3) <- allAssign3$sample

allAssign3$ploidy <-
  as.data.frame(popTable)[rownames(allAssign3), "ploidy"]
allAssign3$species <-
  as.data.frame(popTable)[rownames(allAssign3), "species"]
allAssign3$population <-
  as.data.frame(popTable)[rownames(allAssign3), "population"]
allAssign3$popPloid <- paste(allAssign3$population,
                            allAssign3$ploidy, sep = "")
allAssign3$region <-
  as.data.frame(popTable)[rownames(allAssign3), "area"]

allAssign3$sortA <- paste(allAssign3$ploidy, allAssign3$species,
                         allAssign3$region,
                         allAssign3$population, sep = "") 

allAssign3 <- allAssign3[order(allAssign3$sortA), ]

allAssign3[c("SCCC5", "MOSJ11"), "region" ] <- "tetraploid"

## Four clusters
allAssign4 <- read.csv("structure-all/all4_1.csv", header = FALSE)
allAssign4 <- allAssign4[, -1]
colnames(allAssign4) <- c("sample", "missing", "population",
                         "pop1", "pop2", "pop3", "pop4")
rownames(allAssign4) <- allAssign4$sample

allAssign4$ploidy <-
  as.data.frame(popTable)[rownames(allAssign4), "ploidy"]
allAssign4$species <-
  as.data.frame(popTable)[rownames(allAssign4), "species"]
allAssign4$population <-
  as.data.frame(popTable)[rownames(allAssign4), "population"]
allAssign4$popPloid <- paste(allAssign4$population,
                            allAssign4$ploidy, sep = "")
allAssign4$region <-
  as.data.frame(popTable)[rownames(allAssign4), "area"]

allAssign4$sortA <- paste(allAssign4$ploidy, allAssign4$species,
                         allAssign4$region,
                         allAssign4$population, sep = "") 

allAssign4[c("SCCC5", "MOSJ11"), "region" ] <- "tetraploid"
allAssign4 <- allAssign4[order(allAssign4$sortA), ]


message("all STRUCTURE K")
allK <- read.csv("structure-all/ksummary.csv", header = FALSE)
colnames(allK) <- c("K", "LnProb")

allEvanno <- allK %>% group_by(K) %>%
  summarize(meanEst = mean(LnProb),
            sd = sd(LnProb)) %>%
  mutate(LnP = c(0, diff(meanEst))) %>%
  mutate(LnPP = abs(c(0, diff(LnP)[-1], 0))) %>%
  mutate(deltaK = LnPP/sd)

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

tripPops <- triploids

PopNames(triploids) <- levels(factor(c("South", "North", "C. reticulata")))
PopInfo(triploids)  <- factor(popTable[Samples(triploids), ]$area)

PopNames(tripPops) <- as.character(unique(popTable[Samples(triploids), ]$population))
PopInfo(tripPops)  <- factor(popTable[Samples(tripPops), ]$population)

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

write.Structure(tetraploids, ploidy = 4, file = "structure-all/allPloidy.stru")

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
dipAlleles <- laevAlleles | occiAlleles
northAlleles <- allelefreq["North", -1] > 0  
southAlleles <- allelefreq["South", -1] > 0  
pumiAlleles <- northAlleles | southAlleles
reticAlleles <- allelefreq["C. reticulata", -1] > 0  

lociNames <- substr(colnames(allelefreq), start = 0, stop = 5)


occiPA <- occiAlleles & (! laevAlleles)
laevPA <- laevAlleles & (! occiAlleles)
dipPA <- dipAlleles & !(pumiAlleles)
pumPA <- !dipAlleles & pumiAlleles

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


## Triploid Ho calculation, per population:
## Run the following to get the values in the terminal:

## for(tPop in PopNames(tripPops)) {
##   message(tPop)
##   print(mean(
##     apply(Genotypes(tripPops,
##                     samples = Samples(tripPops, populations = tPop)),
##           1, function(x) sum(sapply(x, length) > 1)/8)))
## }

## for(tPop in PopNames(tripPops)) {
##   message(tPop)
##   print(
##     apply(Genotypes(tripPops,
##                     samples = Samples(tripPops, populations = tPop)),
##           1, function(x) sum(sapply(x, length) > 1)))
## }



## There's a monomorphic marker in the North group, which generates a
## warning here:
options(warn = 1)
simGst <- calcPopDiff(allelefreq, metric = "Gst")
options(warn = 2)
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
sex <- read.table("data/celtis_morph.csv", header = TRUE,
                 sep = "\t", quote = "")

sex[sex == "na"] <- NA

sex <- sex[which(sex$Ploidy == 3),]

sexVars <- c("PetioleL", "LfL", "LfW", "basalveinL", "Lf.widest", "LfTeeth",
            "PedicelL",
            ##"secondaryVeins", 
            "DrupeL")

sex$PetioleL <- as.numeric(sex$PetioleL)
sex$LfL <- as.numeric(sex$LfL)
sex$basalveinL <- as.numeric(sex$basalveinL)
sex$Lf.widest <- as.numeric(sex$Lf.widest)
#sex$secondaryVeins <- as.numeric(sex$secondaryVeins)
sex$DrupeL <- as.numeric(sex$DrupeL)

sex$LW <- sex$LfL/sex$LfW
sex$widestProp <- sex$Lf.widest/sex$LfL
sex$PedPet <- sex$PedicelL/sex$PetioleL
sex$basalVProp <- sex$basalveinL/sex$LfL

sexVars <- c(sexVars, "LW", "widestProp", "PedPet", "basalVProp")

apply(sex[, sexVars], 2, function(x) sum(is.na(x)))
## drop DrupeL, basalveinL, basalVProp: too many missing

sexVars <- sexVars[! sexVars %in% c("DrupeL", "basalveinL", "basalVProp")]

cor(sex[, sexVars], use = "pairwise.complete")
## Dropping highly correlated 'size' variables: LfW, and also
## Lf.widest. Lf.widest not mentioned in manuscript at all.

## Also dropping drupe length: missing in nearly 20% of samples; southern
## might be a little shorter, but lots of overlap, and many fruits are more
## or less distended in pressing, so my confidence in the precise values is
## low.

## and basalVProp, also missing in 13% of samples
sexVars <- c("PetioleL", "LfL", "LfTeeth", "PedicelL", "LW",
            "widestProp", "PedPet") 

sex <- sex[complete.cases(sex[, sexVars]), ]

sex$group <- as.factor(sex$group)
sexLDA <- lda(grouping = sex$group, x = as.data.frame(scale(sex[, sexVars])))

sexCoef <- sexLDA$scaling
row.names(sexCoef) <- c("Petiole Length", "Leaf Length", "Leaf Teeth",
                       "Pedicel Length", "Leaf Length/Width Ratio",
                       "Widest Point of Leaf", "Pedicel/Petiole Ratio")  

sexMeans <- aggregate(sex[, sexVars], by = list(sex$group), mean)
sexSD <- aggregate(sex[, sexVars], by = list(sex$group), sd)

sexCoef <- cbind(sexCoef, t(sexMeans[, -1]), t(sexSD[, -1]))
colnames(sexCoef) <- c("Can. Coef.", "North Mean", "South Mean", "North SD",
                      "South SD")
sexCoef[, 1] <- round(sexCoef[, 1], 2)
sexCoef[, 2:5] <- round(sexCoef[, 2:5], 1)
sexCoef <- data.frame(sexCoef)

sexCoef$North <- sprintf("%.1f ± %.1f", sexCoef[, "North.Mean"],
                        sexCoef[, "North.SD"])
sexCoef$South <- sprintf("%.1f ± %.1f", sexCoef[, "South.Mean"],
                        sexCoef[, "South.SD"])
sexCoef <- sexCoef[, c(1, 6, 7)]

summary(manova(scale(sex[, sexVars]) ~ sex$group), test = "Wilks")

sexPCA <- prcomp(sex[, sexVars], scale. = TRUE, retx = TRUE)
sexPCAEigs <- round(sexPCA$sdev^2 / sum(sexPCA$sdev^2), 3) * 100

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

## Uncomment the following block to regenerate the map data:

## us <- getData("GADM", country = "USA", level = 1, path = "./data/maps/",
##              type = "sf")
## canada <- getData("GADM", country = "CAN", level = 1, path = "./data/maps",
##                  type = "sf")
## mex <- getData("GADM", country = "MEX", level = 1, path = "./data/maps",
##               type = "sf")

## ## update old-style crs:
## st_crs(us) <- st_crs(us)
## st_crs(canada) <- st_crs(canada)
## st_crs(mex) <- st_crs(mex)

## na <- rbind(us, canada, mex)
## na.simp <- st_simplify(na, dTolerance = 0.01)
## laea = CRS("+proj=laea +lat_0=30 +lon_0=-95")
## na.la <- st_transform(na.simp, laea)
na.la.crop <- st_crop(na.la, xmin = -900000, xmax = 2200000,
                      ymin = -500000, ymax = 2600000)
na.la.crop1e3 <- st_simplify(na.la.crop, dTolerance = 1e3)

na.la.crop2 <- st_crop(na.la, xmin = 2e5, xmax = 8e5,
                       ymin = 5e5, ymax = 1.1e6)

na.la.crop2.1e3 <- st_simplify(na.la.crop2, dTolerance = 1e3)

## greatlakes <- st_read("data/maps/greatlakes.shp")
## gl.proj <- st_transform(greatlakes, laea)

## popTableProj <- spTransform(popTable, laea)
save(us, canada, mex, na, na.simp, laea, na.la, greatlakes, gl.proj,
     na.la.crop, popTableProj, file = "data/maps/prepped-maps.Rda")

## After the map data is generated, load it from this file to speed things up:
load("data/maps/prepped-maps.Rda")

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
