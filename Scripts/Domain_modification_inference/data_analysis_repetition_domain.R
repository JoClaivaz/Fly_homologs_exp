'''
Joaquim Claivaz
170520
Script which allows the visualize of informative domain (domain modification not part of domain repetition)
Based on Marta script analysis/plot
'''
setwd('D:/UNIL/Master/Master_Project/Data/domain_architecture_inference')

#library
library(ggplot2) 
library(grid)
#DROME_DROYA
DROME_DROYA <- read.table('all_DROSO/DROME_DROYA/all_ortholog_DROME_DROYA_domain_modif_parsed',header=F, sep="\t")

DROME_DROYA$species <- DROME_DROYA[,2]
colnames(DROME_DROYA) <- c("pair","id","hmm acc", "hmm start", "hmm end", "hmm length", "species")
DROME_DROYA<- aggregate(`hmm acc` ~ pair+species, data = DROME_DROYA, paste, collapse = " ")

DROME_DROYA$status <- with(DROME_DROYA, ave(`hmm acc`, pair,
                                            FUN=function(x) c("Domain loss", "Domain repetitions")[1+(setequal(
                                              unlist(strsplit(x[1],split=' ')),
                                              unlist(strsplit(x[2],split=' '))
                                            )) ]) )
DROME_DROYA_subset<-subset(DROME_DROYA,DROME_DROYA$status=="Domain loss")
DROME_DROYA_domain_loss <-aggregate(species ~ pair, DROME_DROYA_subset, paste, collapse = " ", simplify=FALSE)
DROME_DROYA_domain_loss <- apply(DROME_DROYA_domain_loss,2,as.character)
write.table(DROME_DROYA_domain_loss,file="./all_DROSO/DROME_DROYA/DROME_DROYA_domain_loss",row.names = F,quote = F,col.names = F)

DROME_DROYA$pair_species <- "DROME_DROYA"

keep <- c("pair","status","pair_species")
DROME_DROYA <- DROME_DROYA[keep]
DROME_DROYA <- DROME_DROYA[!duplicated(DROME_DROYA), ]
table(DROME_DROYA$status)


#DROME_DROPS
DROME_DROPS <- read.table('all_DROSO/DROME_DROPS/all_ortholog_DROME_DROPS_domain_modif_parsed',header=F, sep="\t")

DROME_DROPS$species <- DROME_DROPS[,2]
colnames(DROME_DROPS) <- c("pair","id","hmm acc", "hmm start", "hmm end", "hmm length", "species")
DROME_DROPS<- aggregate(`hmm acc` ~ pair+species, data = DROME_DROPS, paste, collapse = " ")

DROME_DROPS$status <- with(DROME_DROPS, ave(`hmm acc`, pair,
                                            FUN=function(x) c("Domain loss", "Domain repetitions")[1+(setequal(
                                              unlist(strsplit(x[1],split=' ')),
                                              unlist(strsplit(x[2],split=' '))
                                            )) ]) )
DROME_DROPS_subset<-subset(DROME_DROPS,DROME_DROPS$status=="Domain loss")
DROME_DROPS_domain_loss <-aggregate(species ~ pair, DROME_DROPS_subset, paste, collapse = " ", simplify=FALSE)
DROME_DROPS_domain_loss <- apply(DROME_DROPS_domain_loss,2,as.character)
write.table(DROME_DROPS_domain_loss,file="./all_DROSO/DROME_DROPS/DROME_DROPS_domain_loss",row.names = F,quote = F,col.names = F)

DROME_DROPS$pair_species <- "DROME_DROPS"

keep <- c("pair","status","pair_species")
DROME_DROPS <- DROME_DROPS[keep]
DROME_DROPS <- DROME_DROPS[!duplicated(DROME_DROPS), ]
table(DROME_DROPS$status)


#DROME_DROAN
DROME_DROAN <- read.table('all_DROSO/DROME_DROAN/all_ortholog_DROME_DROAN_domain_modif_parsed',header=F, sep="\t")

DROME_DROAN$species <- DROME_DROAN[,2]
colnames(DROME_DROAN) <- c("pair","id","hmm acc", "hmm start", "hmm end", "hmm length", "species")
DROME_DROAN<- aggregate(`hmm acc` ~ pair+species, data = DROME_DROAN, paste, collapse = " ")

DROME_DROAN$status <- with(DROME_DROAN, ave(`hmm acc`, pair,
                                            FUN=function(x) c("Domain loss", "Domain repetitions")[1+(setequal(
                                              unlist(strsplit(x[1],split=' ')),
                                              unlist(strsplit(x[2],split=' '))
                                            )) ]) )

DROME_DROAN_subset <- subset(DROME_DROAN,DROME_DROAN$status=="Domain loss")
DROME_DROAN_domain_loss <- aggregate(species ~ pair, DROME_DROAN_subset, paste, collapse = " ", simplify=FALSE)
DROME_DROAN_domain_loss <- apply(DROME_DROAN_domain_loss,2,as.character)
write.table(DROME_DROAN_domain_loss,file="./all_DROSO/DROME_DROAN/DROME_DROAN_domain_loss",row.names = F,quote = F,col.names = F)

DROME_DROAN$pair_species <- "DROME_DROAN"

keep <- c("pair","status","pair_species")
DROME_DROAN <- DROME_DROAN[keep]
DROME_DROAN <- DROME_DROAN[!duplicated(DROME_DROAN), ]
table(DROME_DROAN$status)

#DROME_DROMO
DROME_DROMO <- read.table('all_DROSO/DROME_DROMO/all_ortholog_DROME_DROMO_domain_modif_parsed',header=F, sep="\t")

DROME_DROMO$species <- DROME_DROMO[,2]
colnames(DROME_DROMO) <- c("pair","id","hmm acc", "hmm start", "hmm end", "hmm length", "species")
DROME_DROMO<- aggregate(`hmm acc` ~ pair+species, data = DROME_DROMO, paste, collapse = " ")

DROME_DROMO$status <- with(DROME_DROMO, ave(`hmm acc`, pair,
                                            FUN=function(x) c("Domain loss", "Domain repetitions")[1+(setequal(
                                              unlist(strsplit(x[1],split=' ')),
                                              unlist(strsplit(x[2],split=' '))
                                            )) ]) )

DROME_DROMO_subset <- subset(DROME_DROMO,DROME_DROMO$status=="Domain loss")
DROME_DROMO_domain_loss <- aggregate(species ~ pair, DROME_DROMO_subset, paste, collapse = " ", simplify=FALSE)
DROME_DROMO_domain_loss <- apply(DROME_DROMO_domain_loss,2,as.character)
write.table(DROME_DROMO_domain_loss,file="./all_DROSO/DROME_DROMO/DROME_DROMO_domain_loss",row.names = F,quote = F,col.names = F)

DROME_DROMO$pair_species <- "DROME_DROMO"

keep <- c("pair","status","pair_species")
DROME_DROMO <- DROME_DROMO[keep]
DROME_DROMO <- DROME_DROMO[!duplicated(DROME_DROMO), ]
table(DROME_DROMO$status)

#DROME_DROSI
DROME_DROSI <- read.table('all_DROSO/DROME_DROSI/all_ortholog_DROME_DROSI_domain_modif_parsed',header=F, sep="\t")

DROME_DROSI$species <- DROME_DROSI[,2]
colnames(DROME_DROSI) <- c("pair","id","hmm acc", "hmm start", "hmm end", "hmm length", "species")
DROME_DROSI<- aggregate(`hmm acc` ~ pair+species, data = DROME_DROSI, paste, collapse = " ")

DROME_DROSI$status <- with(DROME_DROSI, ave(`hmm acc`, pair,
                                            FUN=function(x) c("Domain loss", "Domain repetitions")[1+(setequal(
                                              unlist(strsplit(x[1],split=' ')),
                                              unlist(strsplit(x[2],split=' '))
                                            )) ]) )

DROME_DROSI_subset <- subset(DROME_DROSI,DROME_DROSI$status=="Domain loss")
DROME_DROSI_domain_loss <- aggregate(species ~ pair, DROME_DROSI_subset, paste, collapse = " ", simplify=FALSE)
DROME_DROSI_domain_loss <- apply(DROME_DROSI_domain_loss,2,as.character)
write.table(DROME_DROSI_domain_loss,file="./all_DROSO/DROME_DROSI/DROME_DROSI_domain_loss",row.names = F,quote = F,col.names = F)

DROME_DROSI$pair_species <- "DROME_DROSI"

keep <- c("pair","status","pair_species")
DROME_DROSI <- DROME_DROSI[keep]
DROME_DROSI <- DROME_DROSI[!duplicated(DROME_DROSI), ]
table(DROME_DROSI$status)

#DROME_DROVI
DROME_DROVI <- read.table('all_DROSO/DROME_DROVI/all_ortholog_DROME_DROVI_domain_modif_parsed',header=F, sep="\t")

DROME_DROVI$species <- DROME_DROVI[,2]
colnames(DROME_DROVI) <- c("pair","id","hmm acc", "hmm start", "hmm end", "hmm length", "species")
DROME_DROVI<- aggregate(`hmm acc` ~ pair+species, data = DROME_DROVI, paste, collapse = " ")

DROME_DROVI$status <- with(DROME_DROVI, ave(`hmm acc`, pair,
                                            FUN=function(x) c("Domain loss", "Domain repetitions")[1+(setequal(
                                              unlist(strsplit(x[1],split=' ')),
                                              unlist(strsplit(x[2],split=' '))
                                            )) ]) )

DROME_DROVI_subset <- subset(DROME_DROVI,DROME_DROVI$status=="Domain loss")
DROME_DROVI_domain_loss <- aggregate(species ~ pair, DROME_DROVI_subset, paste, collapse = " ", simplify=FALSE)
DROME_DROVI_domain_loss <- apply(DROME_DROVI_domain_loss,2,as.character)
write.table(DROME_DROVI_domain_loss,file="./all_DROSO/DROME_DROVI/DROME_DROVI_domain_loss",row.names = F,quote = F,col.names = F)

DROME_DROVI$pair_species <- "DROME_DROVI"

keep <- c("pair","status","pair_species")
DROME_DROVI <- DROME_DROVI[keep]
DROME_DROVI <- DROME_DROVI[!duplicated(DROME_DROVI), ]
table(DROME_DROVI$status)
