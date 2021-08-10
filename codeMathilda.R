##Test drm + medrm avec mes données

## https://github.com/DoseResponse/medrc
## J'ai téléchargé et installé 'drc' via commandes R 
## J'ai téléchargé et installé 'metafor' via commandes R 
## J'ai téléchargé les fichiers du package : medrc
## J'ai téléchargé les fichiers data de drc
## J'ai installé à la main les 2 (en ayant préalablement installé 'metafor')

library(drc)
library(metafor)
library(medrc)


###Dose-resp models parameters : b : slope, c : lower limit, d : upper limit, e : ED50

################## Tests avec mes données
library("readxl")
setwd("P:/EPI/donnees chaudes/insectarium/Mathilda/Evol Exp/Biotests/sur_G10")
biotestsG10plq <- read_excel("Biotests_plaques_G10.xlsx", "Feuil1")
#On enlève les cases vides
biotestsG10plq<-biotestsG10plq[complete.cases(biotestsG10plq[]),]


########### Comparaison des gradients au sein d'une lignée ########### 

########## Spinosad ########## 

biotestsG10plqSpi<-biotestsG10plq[biotestsG10plq$lignee=="SPI",]
biotestsG10plqSV<-biotestsG10plq[biotestsG10plq$lignee=="SV",]

biotestsG10plqSpiSPI<-biotestsG10plqSpi[biotestsG10plqSpi$insecticide=="SPI" | biotestsG10plqSpi$insecticide=="T",]
biotestsG10plqSVSPI<-biotestsG10plqSV[biotestsG10plqSV$insecticide=="SPI" | biotestsG10plqSV$insecticide=="T",]

biotestsG10plqSpiSPI_0<-biotestsG10plqSpiSPI[biotestsG10plqSpiSPI$LD==0,]
biotestsG10plqSpiSPI_25<-biotestsG10plqSpiSPI[biotestsG10plqSpiSPI$LD==25,]
biotestsG10plqSpiSPI_75<-biotestsG10plqSpiSPI[biotestsG10plqSpiSPI$LD==75,]


#### Pour La LD 0 : comparaison des réplicats ####
modBTSpi_0 <-drm(nb_morts/nb_installes~concentration, 
                 weights=nb_installes, 
                 curveid = replica,
                 data=biotestsG10plqSpiSPI_0,
                 fct=LN.3u(),
                 type="binomial")
modBTSpi_0$replica<-as.factor(modBTSpi_0$replica)

plot(modBTSpi_0, type = "all",
     xlab = "Concentration (ppm)",
     ylab = "Proportion dead",
     main = "Spinosad LD 0",
     legendPos = c(0.2, 1), col = c(1,2,3))

summary(modBTSpi_0)
confint(modBTSpi_0)
compParm(modBTSpi_0, "e", "-")

#### Pour La LD 25 : comparaison des réplicats ####
modBTSpi_25 <-drm(nb_morts/nb_installes~concentration, 
                  weights=nb_installes, 
                  curveid = replica,
                  data=biotestsG10plqSpiSPI_25,
                  fct=LN.3u(),
                  type="binomial")

plot(modBTSpi_25, type = "all",
     xlab = "Concentration (ppm)",
     ylab = "Proportion dead",
     main = "Spinosad LD 25",
     legendPos = c(0.2, 1), col = c(1,2,3))

summary(modBTSpi_25)
confint(modBTSpi_25)
compParm(modBTSpi_25, "e", "-")

#### Pour La LD 75 : comparaison des réplicats ####
modBTSpi_75 <-drm(nb_morts/nb_installes~concentration, 
                  weights=nb_installes, 
                  curveid = replica,
                  data=biotestsG10plqSpiSPI_75,
                  fct=LN.3u(),
                  type="binomial")

plot(modBTSpi_75, type = "all",
     xlab = "Concentration (ppm)",
     ylab = "Proportion dead",
     main = "Spinosad LD 75",
     legendPos = c(0.2, 1), col = c(1,2,3))

summary(modBTSpi_75)
confint(modBTSpi_75)
compParm(modBTSpi_75, "e", "-")

########## Deltaméthrine ########## 
biotestsG10plqDel<-biotestsG10plq[biotestsG10plq$lignee=="DEL",]

biotestsG10plqDelDEL<-biotestsG10plqDel[biotestsG10plqDel$insecticide=="DEL" | biotestsG10plqDel$insecticide=="T",]
biotestsG10plqSVDEL<-biotestsG10plqSV[biotestsG10plqSV$insecticide=="DEL" | biotestsG10plqSV$insecticide=="T",]

biotestsG10plqDelDEL_0<-biotestsG10plqDelDEL[biotestsG10plqDelDEL$LD==0,]
biotestsG10plqDelDEL_50<-biotestsG10plqDelDEL[biotestsG10plqDelDEL$LD==50,]
biotestsG10plqDelDEL_75<-biotestsG10plqDelDEL[biotestsG10plqDelDEL$LD==75,]




#### Pour La LD 0 : comparaison des réplicats ####
modBTDel_0 <-drm(nb_morts/nb_installes~concentration, 
                 weights=nb_installes, 
                 curveid = replica,
                 data=biotestsG10plqDelDEL_0,
                 fct=LN.3u(),
                 type="binomial")

plot(modBTDel_0, type = "all",
     xlab = "Concentration (ppm)",
     ylab = "Proportion dead",
     main = "Deltaméthrine LD 0",
     legendPos = c(0.0005, 1), col = c(1,2,3))

summary(modBTDel_0)
confint(modBTDel_0)
compParm(modBTDel_0, "e", "-")

#### Pour La LD 25 : comparaison des réplicats ####
modBTDel_50 <-drm(nb_morts/nb_installes~concentration, 
                  weights=nb_installes, 
                  curveid = replica,
                  data=biotestsG10plqDelDEL_50,
                  fct=LN.3u(),
                  type="binomial")

plot(modBTDel_50, type = "all",
     xlab = "Concentration (ppm)",
     ylab = "Proportion dead",
     main = "Deltaméthrine LD 50",
     legendPos = c(0.0005, 1), col = c(1,2,3))

summary(modBTDel_50)
confint(modBTDel_50)
compParm(modBTDel_50, "e", "-")

#### Pour La LD 75 : comparaison des réplicats ####
#: erreur car pas assez de points en R1-2-3
#remove R1-2-3 
biotestsG10plqDelDEL_75_crop<-biotestsG10plqDelDEL_75[!biotestsG10plqDelDEL_75$replica=='R1-2-3',]

modBTDel_75 <-drm(nb_morts/nb_installes~concentration, 
                  weights=nb_installes, 
                  curveid = replica,
                  data=biotestsG10plqDelDEL_75_crop,
                  fct=LN.3u(),
                  type="binomial")

plot(modBTDel_75, type = "all",
     xlab = "Concentration (ppm)",
     ylab = "Proportion dead",
     main = "Deltaméthrine LD 75",
     legendPos = c(0.005, 1), col = c(1,2,3))

summary(modBTSpi_75)
confint(modBTSpi_75)
compParm(modBTSpi_75, "e", "-")

########## Rynaxypyr ########## 
biotestsG10plqRy<-biotestsG10plq[biotestsG10plq$lignee=="RY",]

biotestsG10plqRyRY<-biotestsG10plqRy[biotestsG10plqRy$insecticide=="RY" | biotestsG10plqRy$insecticide=="T",]
biotestsG10plqSVRY<-biotestsG10plqSV[biotestsG10plqSV$insecticide=="RY" | biotestsG10plqSV$insecticide=="T",]

biotestsG10plqRyRY_0<-biotestsG10plqRyRY[biotestsG10plqRyRY$LD==0,]
biotestsG10plqRyRY_25<-biotestsG10plqRyRY[biotestsG10plqRyRY$LD==25,]
biotestsG10plqRyRY_75<-biotestsG10plqRyRY[biotestsG10plqRyRY$LD==75,]

#### Pour La LD 0 : comparaison des réplicats ####
modBTRy_0 <-drm(nb_morts/nb_installes~concentration, 
                weights=nb_installes, 
                curveid = replica,
                data=biotestsG10plqRyRY_0,
                fct=LN.3u(),
                type="binomial")

plot(modBTRy_0, type = "all",
     xlab = "Concentration (ppm)",
     ylab = "Proportion dead",
     main = "Rynaxypyr LD 0",
     legendPos = c(0.025, 1), col = c(1,2,3))

summary(modBTRy_0)
confint(modBTRy_0)
compParm(modBTRy_0, "e", "-")


# Pb avec R2 ? 
# --> Regrouper les 3 biotests en un seul ? 
biotestsG10plqRyRY_0_mod<-biotestsG10plqRyRY_0
biotestsG10plqRyRY_0_mod[1, "nb_installes"] = 118
biotestsG10plqRyRY_0_mod[1, "nb_morts"] = 12
biotestsG10plqRyRY_0_mod<- biotestsG10plqRyRY_0_mod[-c(14,21,29,30),] 

modBTRy_0_MOD <-drm(nb_morts/nb_installes~concentration, 
                    weights=nb_installes, 
                    curveid = replica,
                    data=biotestsG10plqRyRY_0_mod,
                    fct=LN.3u(),
                    type="binomial")

plot(modBTRy_0_MOD, type = "all",
     xlab = "Concentration (ppm)",
     ylab = "Proportion dead",
     legendPos = c(0.025, 1), col = c(1,2,3))

summary(modBTRy_0_MOD)
confint(modBTRy_0_MOD)
compParm(modBTRy_0_MOD, "e", "-")


#### Pour La LD 25 : comparaison des réplicats ####
modBTRy_25 <-drm(nb_morts/nb_installes~concentration, 
                 weights=nb_installes, 
                 curveid = replica,
                 data=biotestsG10plqRyRY_25,
                 fct=LN.3u(),
                 type="binomial")

plot(modBTRy_25, type = "all",
     xlab = "Concentration (ppm)",
     ylab = "Proportion dead",
     main = "Rynaxypyr LD 25",
     legendPos = c(0.025, 1), col = c(1,2,3))

summary(modBTRy_25)
confint(modBTRy_25)
compParm(modBTRy_25, "e", "-")

#### Pour La LD 75 : comparaison des réplicats ####
modBTRy_75 <-drm(nb_morts/nb_installes~concentration, 
                 weights=nb_installes, 
                 curveid = replica,
                 data=biotestsG10plqRyRY_75,
                 fct=LN.3u(),
                 type="binomial")

plot(modBTRy_75, type = "all",
     xlab = "Concentration (ppm)",
     ylab = "Proportion dead",
     main = "Rynaxypyr LD 75",
     legendPos = c(0.025, 1), col = c(1,2,3))

summary(modBTSpi_75)
confint(modBTSpi_75)
compParm(modBTSpi_75, "e", "-")


########### 


########### Comparaison entre doses ########### 
########## Spinosad ########## 

## drm seul ##
#### Correspond à l'exemple fluoranthene dans le book dose response (binomial + multiple DR curves)
modBTSpi<-drm(nb_morts/nb_installes~concentration, 
              weights=nb_installes, 
              curveid = LD,
              data=biotestsG10plqSpiSPI,
              fct=LN.3u(),
              type="binomial")

summary(modBTSpi)
confint(modBTSpi)
compParm(modBTSpi, "e", "-")
EDcomp(modBTSpi, c(50,50))

plot(modBTSpi, type = "all",
     xlab = "Concentration (ppm)",
     ylab = "Proportion dead",
     main = "Spinosad LD 0 vs 25 vs 75",
     legendPos = c(0.2, 1), col = c(1,2,3))


## Ajout de la SV
biotestsG10plqSpiSPISV<-merge(biotestsG10plqSpiSPI,biotestsG10plqSVSPI,all=T)
modBTSpiSv<-drm(nb_morts/nb_installes~concentration, 
                weights=nb_installes, 
                curveid = LD,
                data=biotestsG10plqSpiSPISV,
                fct=LN.3u(),
                type="binomial")

summary(modBTSpiSv)
confint(modBTSpiSv)
compParm(modBTSpiSv, "e", "-")
EDcomp(modBTSpiSv, c(50,50))

plot(modBTSpiSv, type = "all",
     xlab = "Concentration (ppm)",
     ylab = "Proportion dead",
     main = "Spinosad LD 0 vs 25 vs 75 vs SV",
     legendPos = c(0.2, 1), col = c(1,2,3))


########## Deltaméthrine ########## 

## drm seul ##
# Correspond à l'exemple fluoranthene dans le book dose response (binomial + multiple DR curves)
biotestsG10plqDelDEL_crop<-biotestsG10plqDelDEL[!biotestsG10plqDelDEL$replica=='R1-2-3',]

modBTDel<-drm(nb_morts/nb_installes~concentration, 
              weights=nb_installes, 
              curveid = LD,
              data=biotestsG10plqDelDEL_crop,
              fct=LN.3u(),
              type="binomial")

summary(modBTDel)
confint(modBTDel)
compParm(modBTDel, "e", "-")
EDcomp(modBTDel, c(50,50))


plot(modBTDel, type = "all",
     xlab = "Concentration (ppm)",
     ylab = "Proportion dead",
     main = "Deltaméthrine LD 0 vs 50 vs 75",
     legendPos = c(0.0005, 1), col = c(1,2,3))


## Ajout de la SV
biotestsG10plqDelDELSV<-merge(biotestsG10plqDelDEL,biotestsG10plqSVDEL,all=T)
modBTDelSv<-drm(nb_morts/nb_installes~concentration, 
                weights=nb_installes, 
                curveid = LD,
                data=biotestsG10plqDelDELSV,
                fct=LN.3u(),
                type="binomial")

summary(modBTDelSv)
confint(modBTDelSv)
compParm(modBTDelSv, "e", "-")
EDcomp(modBTDelSv, c(50,50))

plot(modBTDelSv, type = "all",
     xlab = "Concentration (ppm)",
     ylab = "Proportion dead",
     main = "Deltaméthrine LD 0 vs 50 vs 75 vs SV",
     legendPos = c(0.2, 1), col = c(1,2,3))


########## Rynaxypyr ########## 

## drm seul ##
# Correspond à l'exemple fluoranthene dans le book dose response
modBTRy<-drm(nb_morts/nb_installes~concentration, 
             weights=nb_installes, 
             curveid = LD,
             data=biotestsG10plqRyRY,
             fct=LN.3u(),
             type="binomial")


summary(modBTRy)
confint(modBTRy)
compParm(modBTRy, "e", "-")
EDcomp(modBTRy, c(50,50))


plot(modBTRy, type = "all",
     xlab = "Concentration (ppm)",
     ylab = "Proportion dead",
     main = "Rynaxypyr LD 0 vs 25 vs 75",
     legendPos = c(0.025, 1), col = c(1,2,3))


## Ajout de la SV
biotestsG10plqRyRYSV<-merge(biotestsG10plqRyRY,biotestsG10plqSVRY,all=T)
modBTRySv<-drm(nb_morts/nb_installes~concentration, 
               weights=nb_installes, 
               curveid = LD,
               data=biotestsG10plqRyRYSV,
               fct=LN.3u(),
               type="binomial")

summary(modBTRySv)
confint(modBTRySv)
compParm(modBTRySv, "e", "-")
EDcomp(modBTRySv, c(50,50))

plot(modBTRySv, type = "all",
     xlab = "Concentration (ppm)",
     ylab = "Proportion dead",
     main = "Rynaxypyr LD 0 vs 25 vs 75 vs SV",
     legendPos = c(0.025, 1), col = c(1,2,3))





########### 

#### Tests medrm ~...~ ####

##medrm SPINO
# LD 0
metamodSpiSPI_0 <- metadrm(nb_morts/nb_installes~concentration,
                           data = biotestsG10plqSpiSPI_0,
                           fct = LN.3u(),
                           ind = replica,
                           struct = "UN")
summary(metamodSpiSPI_0)


# LD 25
metamodSpiSPI_25 <- metadrm(nb_morts/nb_installes~concentration,
                            data = biotestsG10plqSpiSPI_25,
                            fct = LN.3u(),
                            ind = replica,
                            struct = "UN")
summary(metamodSpiSPI_25)

# LD 75
metamodSpiSPI_75 <- metadrm(nb_morts/nb_installes~concentration,
                            data = biotestsG10plqSpiSPI_75,
                            fct = LN.3u(),
                            ind = replica,
                            struct = "UN")
summary(metamodSpiSPI_75)

## Toutes les LD avec rep en effet aléatoire + distinction des groupes avec LD

biotestsG10plqSpiSPI$LD<-as.factor(biotestsG10plqSpiSPI$LD)
biotestsG10plqSpiSPI$replica<-as.factor(biotestsG10plqSpiSPI$replica)


#metamodSpiSPI <- metadrm(nb_morts/nb_installes~concentration,
#                         data = biotestsG10plqSpiSPI,
#                         fct = LN.3u(),
#                         ind = replica,
#                         cid2 = LD, 
#                         struct = "UN")
#summary(metamodSpiSPI)
##erreur ? j'utilise curveid à la place ?

metamodSpiSPI <- metadrm(nb_morts/nb_installes~concentration,
                         data = biotestsG10plqSpiSPI,
                         fct = LL.3u(),
                         ind = replica,
                         curveid = LD, 
                         struct = "UN")
summary(metamodSpiSPI)

EDcomp(metamodSpiSPI,
       percVec = c(50),
       percMat = rbind(c(1, 1)),
       interval = "delta")

modBTspi2 <- medrm(nb_morts/nb_installes~concentration,
                   curveid = b + c + e ~ LD,
                   data = biotestsG10plqSpiSPI,
                   fct = LN.3u(),
                   random = b + c + e ~ 1|replica)

