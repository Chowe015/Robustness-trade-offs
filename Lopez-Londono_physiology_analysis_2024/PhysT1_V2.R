#Robustness trade-offs: an applied ecological approach to understand niche partitioning 
#in sibling species of reef corals

#Lopez-Londoño and Howe, et al.

# set working directory to  folder containing required excel files
getwd()

library("xlsx")
library("lsmeans")

# Import Data
PhysT1 <- read.xlsx("AllPhysUSVI_T1.xlsx", sheetName="AllPhysUSVI_T1")
PhysMeans <- read.xlsx("AllPhysUSVI_T1.xlsx", sheetName="Means")
FvFm <- read.xlsx("FvFmUSVI_T1.xlsx", sheetName="FvFmT1")
dFFmSep <- read.xlsx("FvFmUSVI_T1.xlsx", sheetName="dFFm")
FvFmSep <- read.xlsx("FvFmUSVI_T1.xlsx", sheetName="FvFm")
QmSep <- read.xlsx("FvFmUSVI_T1.xlsx", sheetName="Qm")

FvFmGrad <- read.xlsx("PSIIgradient.xlsx", sheetName="FvFm")
dFFmGrad <- read.xlsx("PSIIgradient.xlsx", sheetName="dFFm.")
QmGrad <- read.xlsx("PSIIgradient.xlsx", sheetName="Qm")

temp <- read.xlsx("TempLongFormat.xlsx", sheetName="LongFormat")
irrad <- read.xlsx("Both sites-8 days.xlsx", sheetName="DLI")

#Test for normality
set.seed(666)
norm.Chla <- shapiro.test(PhysT1$Chlam2)        #Normal distribution (W=0.96; p=0.1577)
norm.SymbDens <- shapiro.test(PhysT1$SymbDens)  #Normal distribution (W=0.97, p=0.274)
norm.Ci <- shapiro.test(PhysT1$Ci)              #Normal distribution (W=0.98; p=0.7454)
norm.Prot <- shapiro.test(PhysT1$Prot)          #Normal distribution (W=0.98; p=0.7536)
norm.APAR <- shapiro.test(PhysT1$A.PAR)         #NON-normal distribution (W=0.91; p=0.0012)
norm.logAPAR <- shapiro.test(log(PhysT1$A.PAR)) #NON-normal distribution
norm.sqrtAPAR <- shapiro.test(sqrt(PhysT1$A.PAR))   #NON-normal distribution
norm.a <- shapiro.test(PhysT1$a.)               #NON-normal distribution (W=0.94; p=0.0265)
norm.aSym <- shapiro.test(PhysT1$a.Symb)        #Normal distribution (W=0.95; p=0.0608)
Norm.Alfa <- shapiro.test(PhysT1$Alfa)          #Normal distribution (W=0.95; p=0.0777)
norm.Ec <- shapiro.test(PhysT1$Ec)              #Normal distribution (W=0.95; p=0.0539)
norm.Ec <- shapiro.test(PhysT1$Ec)              #Normal distribution (W=0.95; p=0.0724)
norm.Ek <- shapiro.test(PhysT1$Ek)              #Normal distribution (W=0.97; p=0.2894)
norm.Rd <- shapiro.test(PhysT1$RdAVG)           #Normal distribution (W=0.97; p=0.3057)
norm.Pmax <- shapiro.test(PhysT1$Pmaxg)         #Normal distribution (W=0.97; p=0.4392)
norm.PR <- shapiro.test(PhysT1$P.R)             #Normal distribution (W=0.96; p=0.1236)
norm.MQY <- shapiro.test(PhysT1$MQY)            #Normal distribution (W=0.97; p=0.3116)
norm.MRQ <- shapiro.test(PhysT1$MRQ)            #Normal distribution (W=0.96; p=0.1368)
norm.FvFm <- shapiro.test(FvFm$Fv.Fm)           #Normal distribution (W=0.98; p=0.0586)
norm.dFFm <- shapiro.test(FvFm$dF.Fm.)          #NON-normal distribution (W=0.94; p<0.001)
norm.logdFFm <- shapiro.test(log(FvFm$dF.Fm.))  #NON-normal distribution
norm.sqrtdFFm <- shapiro.test(sqrt(FvFm$dF.Fm.))   #NON-normal distribution
norm.Qm <- shapiro.test(FvFm$Qm)                #NON-normal distribution (W=0.93; p<0.001)
norm.FvFmGrad <- shapiro.test(FvFmGrad$Yield)   #Normal distribution (W=0.99, p=0.9386)
norm.dFFmGrad <- shapiro.test(dFFmGrad$Yield)   #NON-normal distribution (W=0.91, p=0.00085)
norm.QmGrad <- shapiro.test(QmGrad$Yield)       #NON-normal distribution (W=0.90, p=0.00028)
norm.Temp <- shapiro.test(temp$Temp)            #NON-normal distribution (W=0.97, p=3.71e-08)
norm.Irrad <- shapiro.test(irrad$DLI)           #NON-normal distribution (W=0.84, p=0.01105)
norm.IrradMx <- shapiro.test(irrad$MaxE)        #NON-normal distribution (W=0.84, p=0.004976)


#MORTALITY
#Fisher's exact test for mortality
survivorship <- data.frame(
    "Shallow" = c(80.85, 80.95),
    "Deep" = c(74.47, 76.19),
    row.names = c("Ofav", "Ofra"),
    stringsAsFactors = FALSE
)
colnames(survivorship) <- c("Shallow", "Deep")
survivorship
chisq.test(survivorship)$expected
fishr <- fisher.test(survivorship)
fishr
#p-value=0.9103, differences in mortality were non-significant


##STATISTICAL ANALYSIS FOR PARAMETRIC DATA##

#CHLOROPHYLL a FLUORESCENCE#
#Fv/Fm
FvFm <- read.xlsx("FvFmUSVI_T1.xlsx", sheetName="FvFmT1")
FvFm.Sep.SH <- t.test(FvFmSep$Ofav.SH, FvFmSep$Ofra.SH, paired=FALSE, exact=FALSE)
FvFm.Sep.SH #FvFm is NOT significantly different among species in the shallow site, p=0.3186
FvFm.Sep.DP <- t.test(FvFmSep$Ofav.DP, FvFmSep$Ofra.DP, paired=FALSE, exact=FALSE)
FvFm.Sep.DP #FvFm is significantly different among species in the deep site, p=6.664e-05***
FvFm.Sep.Ofav <- t.test(FvFmSep$Ofav.SH, FvFmSep$Ofav.DP, paired=FALSE, exact=FALSE)
FvFm.Sep.Ofav #FvFm is significantly different among depths for Ofav, p=0.01515*
FvFm.Sep.Ofra <- t.test(FvFmSep$Ofra.SH, FvFmSep$Ofra.DP, paired=FALSE, exact=FALSE)
FvFm.Sep.Ofra #FvFm is significantly different among depths for Ofra, p=4.98e-14***


#STRUCTURAL TRAITS#
#Chlorophyll a density
Chla <- read.xlsx("AllPhysUSVI_T1.xlsx", sheetName="Chla")
View(Chla)
Chla.Sep.SH <- t.test(Chla$Ofav.SH, Chla$Ofra.SH, paired=FALSE, exact=FALSE)
Chla.Sep.SH #Chla is significantly different among species in the shallow site, p=0.0006657***
Chla.Sep.DP <- t.test(Chla$Ofav.DP, Chla$Ofra.DP, paired=FALSE, exact=FALSE)
Chla.Sep.DP #Chla is significantly different among species in the deep site, p=0.0005518***
Chla.Sep.Ofav <- t.test(Chla$Ofav.SH, Chla$Ofav.DP, paired=FALSE, exact=FALSE)
Chla.Sep.Ofav #Chla is significantly different among depths for Ofav, p=0.004978**
Chla.Sep.Ofra <- t.test(Chla$Ofra.SH, Chla$Ofra.DP, paired=FALSE, exact=FALSE)
Chla.Sep.Ofra #Chla is NOT significantly different among depths for Ofra, p=0.1139

#Symbiont density
SymbDen <- read.xlsx("AllPhysUSVI_T1.xlsx", sheetName="SymbDen")
View(SymbDen)
SymDen.Sep.SH <- t.test(SymbDen$Ofav.SH, SymbDen$Ofra.SH, paired=FALSE, exact=FALSE)
SymDen.Sep.SH #SymbDen is significantly different among species in the shallow site, p=0.002811**
SymDen.Sep.DP <- t.test(SymbDen$Ofav.DP, SymbDen$Ofra.DP, paired=FALSE, exact=FALSE)
SymDen.Sep.DP #Symb.Dens. is significantly different among species in the deep site, p=0.003495**
SymDen.Sep.Ofav <- t.test(SymbDen$Ofav.SH, SymbDen$Ofav.DP, paired=FALSE, exact=FALSE)
SymDen.Sep.Ofav #Symb.Dens. is NOT significantly different among depths for Ofav, p=0.943
SymDen.Sep.Ofra <- t.test(SymbDen$Ofra.SH, SymbDen$Ofra.DP, paired=FALSE, exact=FALSE)
SymDen.Sep.Ofra #Symb.Dens. is NOT significantly different among depths for Ofra, p=0.5153

#Chlorophyll a per symbiont, Ci
Ci <- read.xlsx("AllPhysUSVI_T1.xlsx", sheetName="Ci")
Ci.Sep.SH <- t.test(Ci$Ofav.SH, Ci$Ofra.SH, paired=FALSE, exact=FALSE); Ci.Sep.SH
Ci.Sep.SH #Ci is NOT significantly different among species in the shallow site, p=0.06438
Ci.Sep.DP <- t.test(Ci$Ofav.DP, Ci$Ofra.DP, paired=FALSE, exact=FALSE); Ci.Sep.DP
Ci.Sep.DP #Ci is significantly different among species in the deep site, p=0.000622***
Ci.Sep.Ofav <- t.test(Ci$Ofav.SH, Ci$Ofav.DP, paired=FALSE, exact=FALSE); Ci.Sep.Ofav
Ci.Sep.Ofav #Ci is significantly different among depths for Ofav, p=0.0001936***
Ci.Sep.Ofra <- t.test(Ci$Ofra.SH, Ci$Ofra.DP, paired=FALSE, exact=FALSE); Ci.Sep.Ofra
Ci.Sep.Ofra #Ci is significantly different among depths for Ofra, p=0.01065*

#Total host tissue protein
Prot <- read.xlsx("AllPhysUSVI_T1.xlsx", sheetName="Prot")
Prot.Sep.SH <- t.test(Prot$Ofav.SH, Prot$Ofra.SH, paired=FALSE, exact=FALSE); Prot.Sep.SH
Prot.Sep.SH #Prot is NOT significantly different among species in the shallow site, p=0.3717
Prot.Sep.DP <- t.test(Prot$Ofav.DP, Prot$Ofra.DP, paired=FALSE, exact=FALSE); Prot.Sep.DP
Prot.Sep.DP #Prot is NOT significantly different among species in the deep site, p=0.1977
Prot.Sep.Ofav <- t.test(Prot$Ofav.SH, Prot$Ofav.DP, paired=FALSE, exact=FALSE); Prot.Sep.Ofav
Prot.Sep.Ofav #Prot is NOT significantly different among depths for Ofav, p=0.6605
Prot.Sep.Ofra <- t.test(Prot$Ofra.SH, Prot$Ofra.DP, paired=FALSE, exact=FALSE); Prot.Sep.Ofra
Prot.Sep.Ofra #Prot is NOT significantly different among depths for Ofra, p=0.3667


#OPTICAL PROPERTIES#
#Symbiont specific absorption coefficient 
a.Sym <- read.xlsx("AllPhysUSVI_T1.xlsx", sheetName="a.Sym")
a.Sym.Sep.SH <- t.test(a.Sym$Ofav.SH, a.Sym$Ofra.SH, paired=FALSE, exact=FALSE); a.Sym.Sep.SH
a.Sym.Sep.SH #a.Sym is significantly different among species in the shallow site, p=0.01264*
a.Sym.Sep.DP <- t.test(a.Sym$Ofav.DP, a.Sym$Ofra.DP, paired=FALSE, exact=FALSE); a.Sym.Sep.DP
a.Sym.Sep.DP #a.Sym is significantly different among species in the deep site, p=0.005562**
a.Sym.Sep.Ofav <- t.test(a.Sym$Ofav.SH, a.Sym$Ofav.DP, paired=FALSE, exact=FALSE); a.Sym.Sep.Ofav
a.Sym.Sep.Ofav #a.Sym is significantly different among depths for Ofav, p=0.01547*
a.Sym.Sep.Ofra <- t.test(a.Sym$Ofra.SH, a.Sym$Ofra.DP, paired=FALSE, exact=FALSE); a.Sym.Sep.Ofra
a.Sym.Sep.Ofra #a.Sym is NOT significantly different among depths for Ofra, p=0.06067


#PHOTOSYNTHETIC PARAMETERS#
#Photosynthetic efficiency, alpha
Alfa <- read.xlsx("AllPhysUSVI_T1.xlsx", sheetName="Alfa")
Alfa.Sep.SH <- t.test(Alfa$Ofav.SH, Alfa$Ofra.SH, paired=FALSE, exact=FALSE); Alfa.Sep.SH
Alfa.Sep.SH #Alfa is NOT significantly different among species in the shallow site, p=0.7808
Alfa.Sep.DP <- t.test(Alfa$Ofav.DP, Alfa$Ofra.DP, paired=FALSE, exact=FALSE); Alfa.Sep.DP
Alfa.Sep.DP #Alfa is NOT significantly different among species in the deep site, p=0.5686
Alfa.Sep.Ofav <- t.test(Alfa$Ofav.SH, Alfa$Ofav.DP, paired=FALSE, exact=FALSE); Alfa.Sep.Ofav
Alfa.Sep.Ofav #Alfa is significantly different among depths for Ofav, p=0.007419**
Alfa.Sep.Ofra <- t.test(Alfa$Ofra.SH, Alfa$Ofra.DP, paired=FALSE, exact=FALSE); Alfa.Sep.Ofra
Alfa.Sep.Ofra #Alfa is significantly different among depths for Ofra, p=0.001396**

#Compensating irradiance, Ec
Ec <- read.xlsx("AllPhysUSVI_T1.xlsx", sheetName="Ec")
Ec.Sep.SH <- t.test(Ec$Ofav.SH, Ec$Ofra.SH, paired=FALSE, exact=FALSE); Ec.Sep.SH
Ec.Sep.SH #Ec is significantly different among species in the shallow site, p=0.002741**
Ec.Sep.DP <- t.test(Ec$Ofav.DP, Ec$Ofra.DP, paired=FALSE, exact=FALSE); Ec.Sep.DP
Ec.Sep.DP #Ec is NOT significantly different among species in the deep site, p=0.3175
Ec.Sep.Ofav <- t.test(Ec$Ofav.SH, Ec$Ofav.DP, paired=FALSE, exact=FALSE); Ec.Sep.Ofav
Ec.Sep.Ofav #Ec is significantly different among depths for Ofav, p=8.696e-06***
Ec.Sep.Ofra <- t.test(Ec$Ofra.SH, Ec$Ofra.DP, paired=FALSE, exact=FALSE); Ec.Sep.Ofra
Ec.Sep.Ofra #Ec is significantly different among depths for Ofra, p=0.0005857***

#Saturating irradiance, Ek
Ek <- read.xlsx("AllPhysUSVI_T1.xlsx", sheetName="Ek")
Ek.Sep.SH <- t.test(Ek$Ofav.SH, Ek$Ofra.SH, paired=FALSE, exact=FALSE); Ek.Sep.SH
Ek.Sep.SH #Ek is significantly different among species in the shallow site, p=0.004152**
Ek.Sep.DP <- t.test(Ek$Ofav.DP, Ek$Ofra.DP, paired=FALSE, exact=FALSE); Ek.Sep.DP
Ek.Sep.DP #Ek is NOT significantly different among species in the deep site, p=0.6635
Ek.Sep.Ofav <- t.test(Ek$Ofav.SH, Ek$Ofav.DP, paired=FALSE, exact=FALSE); Ek.Sep.Ofav
Ek.Sep.Ofav #Ek is significantly different among depths for Ofav, p=1.103e-06***
Ek.Sep.Ofra <- t.test(Ek$Ofra.SH, Ek$Ofra.DP, paired=FALSE, exact=FALSE); Ek.Sep.Ofra
Ek.Sep.Ofra #Ek is significantly different among depths for Ofra, p=3.914e-07***

#Respiration rates per unit area, Rd
Rd <- read.xlsx("AllPhysUSVI_T1.xlsx", sheetName="Rd")
Rd.Sep.SH <- t.test(Rd$Ofav.SH, Rd$Ofra.SH, paired=FALSE, exact=FALSE); Rd.Sep.SH
Rd.Sep.SH #Rd is significantly different among species in the shallow site, p=0.02005*
Rd.Sep.DP <- t.test(Rd$Ofav.DP, Rd$Ofra.DP, paired=FALSE, exact=FALSE); Rd.Sep.DP
Rd.Sep.DP #Rd is NOT significantly different among species in the deep site, p=0.248
Rd.Sep.Ofav <- t.test(Rd$Ofav.SH, Rd$Ofav.DP, paired=FALSE, exact=FALSE); Rd.Sep.Ofav
Rd.Sep.Ofav #Rd is significantly different among depths for Ofav, p=0.00422**
Rd.Sep.Ofra <- t.test(Rd$Ofra.SH, Rd$Ofra.DP, paired=FALSE, exact=FALSE); Rd.Sep.Ofra
Rd.Sep.Ofra #Rd is NOT significantly different among depths for Ofra, p=0.09838

#Respiration rates per symbiont, Rsym
Rd.Sym <- read.xlsx("AllPhysUSVI_T1.xlsx", sheetName="Rd.Sym")
Rd.Sym.Sep.SH <- t.test(Rd.Sym$Ofav.SH, Rd.Sym$Ofra.SH, paired=FALSE, exact=FALSE); Rd.Sym.Sep.SH
Rd.Sym.Sep.SH #Rd per symb is NOT significantly different among species in the shallow site, p=0.2806
Rd.Sym.Sep.DP <- t.test(Rd.Sym$Ofav.DP, Rd.Sym$Ofra.DP, paired=FALSE, exact=FALSE); Rd.Sym.Sep.DP
Rd.Sym.Sep.DP #Rd per symb is NOT significantly different among species in the deep site, p=0.146
Rd.Sym.Sep.Ofav <- t.test(Rd.Sym$Ofav.SH, Rd.Sym$Ofav.DP, paired=FALSE, exact=FALSE); Rd.Sym.Sep.Ofav
Rd.Sym.Sep.Ofav #Rd per symb is significantly different among depths for Ofav, p=0.007913**
Rd.Sym.Sep.Ofra <- t.test(Rd.Sym$Ofra.SH, Rd.Sym$Ofra.DP, paired=FALSE, exact=FALSE); Rd.Sym.Sep.Ofra
Rd.Sym.Sep.Ofra #Rd per symb is NOT significantly different among depths for Ofra, p=0.4566

#Maximum photosynthesis per unit area, Pmax
Pmax <- read.xlsx("AllPhysUSVI_T1.xlsx", sheetName="Pmax")
Pmax.Sep.SH <- t.test(Pmax$Ofav.SH, Pmax$Ofra.SH, paired=FALSE, exact=FALSE); Pmax.Sep.SH
Pmax.Sep.SH #Pmax is significantly different among species in the shallow site, p=0.01773*
Pmax.Sep.DP <- t.test(Pmax$Ofav.DP, Pmax$Ofra.DP, paired=FALSE, exact=FALSE); Pmax.Sep.DP
Pmax.Sep.DP #Pmax is NOT significantly different among species in the deep site, p=0.4959
Pmax.Sep.Ofav <- t.test(Pmax$Ofav.SH, Pmax$Ofav.DP, paired=FALSE, exact=FALSE); Pmax.Sep.Ofav
Pmax.Sep.Ofav #Pmax is significantly different among depths for Ofav, p=0.001376**
Pmax.Sep.Ofra <- t.test(Pmax$Ofra.SH, Pmax$Ofra.DP, paired=FALSE, exact=FALSE); Pmax.Sep.Ofra
Pmax.Sep.Ofra #Pmax is NOT significantly different among depths for Ofra, p=0.08463

#Maximum photosynthesis per symbiont, Psym
Pmax.Sym <- read.xlsx("AllPhysUSVI_T1.xlsx", sheetName="Pmax.Sym")
Pmax.Sym.Sep.SH <- t.test(Pmax.Sym$Ofav.SH, Pmax.Sym$Ofra.SH, paired=FALSE, exact=FALSE); Pmax.Sym.Sep.SH
Pmax.Sym.Sep.SH #Pmax per symb is significantly different among species in the shallow site, p=0.006386**
Pmax.Sym.Sep.DP <- t.test(Pmax.Sym$Ofav.DP, Pmax.Sym$Ofra.DP, paired=FALSE, exact=FALSE); Pmax.Sym.Sep.DP
Pmax.Sym.Sep.DP #Pmax per symb is NOT significantly different among species in the deep site, p=0.3092
Pmax.Sym.Sep.Ofav <- t.test(Pmax.Sym$Ofav.SH, Pmax.Sym$Ofav.DP, paired=FALSE, exact=FALSE); Pmax.Sym.Sep.Ofav
Pmax.Sym.Sep.Ofav #Pmax per symb is significantly different among depths for Ofav, p=0.04524*
Pmax.Sym.Sep.Ofra <- t.test(Pmax.Sym$Ofra.SH, Pmax.Sym$Ofra.DP, paired=FALSE, exact=FALSE); Pmax.Sym.Sep.Ofra
Pmax.Sym.Sep.Ofra #Pmax per symb is NOT significantly different among depths for Ofra, p=0.06681

#Minimum quantum requirement, 1/O
MRQ <- read.xlsx("AllPhysUSVI_T1.xlsx", sheetName="MQR")
MRQ.Sep.SH <- t.test(MRQ$Ofav.SH, MRQ$Ofra.SH, paired=FALSE, exact=FALSE); MRQ.Sep.SH
MRQ.Sep.SH #MRQ is NOT significantly different among species in the shallow site, p=0.2952
MRQ.Sep.DP <- t.test(MRQ$Ofav.DP, MRQ$Ofra.DP, paired=FALSE, exact=FALSE); MRQ.Sep.DP
MRQ.Sep.DP #MRQ is NOT significantly different among species in the deep site, p=0.8131
MRQ.Sep.Ofav <- t.test(MRQ$Ofav.SH, MRQ$Ofav.DP, paired=FALSE, exact=FALSE); MRQ.Sep.Ofav
MRQ.Sep.Ofav #MRQ is significantly different among depths for Ofav, p=0.03259*
MRQ.Sep.Ofra <- t.test(MRQ$Ofra.SH, MRQ$Ofra.DP, paired=FALSE, exact=FALSE); MRQ.Sep.Ofra
MRQ.Sep.Ofra #MRQ is significantly different among depths for Ofra, p=0.02097*


##STATISTICAL ANALYSES FOR NON-PARAMETRIC DATA##

#ENVIRONMENTAL CONDITIONS#
#Temperature
MannW.Temp <- wilcox.test(temp$TempSH, temp$TempDP)
MannW.Temp2 <- wilcox.test(Temp ~ Depth, data=temp)
MannW.Temp2 #Temperature is significantly different among sites, W=25115, p=5.93e-05***

#Irradiance
MannW.Irrad <- wilcox.test(irrad$DLI.SH, irrad$DLI.DP)
MannW.Irrad2 <- wilcox.test(DLI ~ Depth, data=irrad)
MannW.Irrad2 #Irradiance is significanlty different among sites, W=63, p=0.0001748***
MannW.MaxE <- wilcox.test(irrad$Max.SH, irrad$Max.DP)
MannW.MaxE2 <- wilcox.test(MaxE ~ Depth2, data=irrad)
MannW.MaxE2 #Max Irrad is significanlty different among sites, W=80, p=4.57e-05***


#CHLOROPHYLL a FLUORESCENCE#
#Effective quantum yield of PSII at noon, DF/Fm'
MannW.dFFm.Sh <- wilcox.test(dFFmSep$Ofav.SH, dFFmSep$Ofra.SH, paired=FALSE, exact=FALSE)
MannW.dFFm.Sh #dFFm' is significantly different among species in the shallow site, p=3.15e-13***
MannW.dFFm.Dp <- wilcox.test(dFFmSep$Ofav.DP, dFFmSep$Ofra.DP, paired=FALSE, exact=FALSE)
MannW.dFFm.Dp #dFFm' is significantly different among species in the deep site, p=6.29e-07***
MannW.dFFm.Ofav <- wilcox.test(dFFmSep$Ofav.SH, dFFmSep$Ofav.DP, paired=FALSE, exact=FALSE)
MannW.dFFm.Ofav #dFFm' is significantly different among sites for Ofav, p=2.46e-14***
MannW.dFFm.Ofra <- wilcox.test(dFFmSep$Ofra.SH, dFFmSep$Ofra.DP, paired=FALSE, exact=FALSE)
MannW.dFFm.Ofra #dFFm' is significantly different among sites for Ofra, p=2.2e-16***

#Maximum excitation pressure over PSII at noon, Qm
MannW.Qm.Sh <- wilcox.test(QmSep$Ofav.SH, QmSep$Ofra.SH, paired=FALSE, exact=FALSE)
MannW.Qm.Sh #Qm is significantly different among species in the shallow site, p=1.03e-12***
MannW.Qm.Dp <- wilcox.test(QmSep$Ofav.DP, QmSep$Ofra.DP, paired=FALSE, exact=FALSE)
MannW.Qm.Dp #Qm is NOT significantly different among species in the deep site, p=0.0896
MannW.Qm.Ofav <- wilcox.test(QmSep$Ofav.SH, QmSep$Ofav.DP, paired=FALSE, exact=FALSE)
MannW.Qm.Ofav #Qm is significantly different among sites for Ofav, p=4.75e-11***
MannW.Qm.Ofra <- wilcox.test(QmSep$Ofra.SH, QmSep$Ofra.DP, paired=FALSE, exact=FALSE)
MannW.Qm.Ofra #Qm is significantly different among sites for Ofra, p=4.76e-15***


#OPTICAL PROPERTIES#
#Absorptance, A-PAR
APAR.Sep <- read.xlsx("AllPhysUSVI_T1.xlsx", sheetName="APAR")
MannW.APAR.Sh <- wilcox.test(APAR.Sep$Ofav.SH, APAR.Sep$Ofra.SH, paired=FALSE, exact=FALSE)
MannW.APAR.Sh #APAR is NOT significantly different among species in the shallow site, p=1
MannW.APAR.Dp <- wilcox.test(APAR.Sep$Ofav.DP, APAR.Sep$Ofra.DP, paired=FALSE, exact=FALSE)
MannW.APAR.Dp #APAR is NOT significantly different among species in the deep site, p=0.1858
MannW.APAR.Ofav <- wilcox.test(APAR.Sep$Ofav.SH, APAR.Sep$Ofav.DP, paired=FALSE, exact=FALSE)
MannW.APAR.Ofav #APAR is significantly different among sites for Ofav, p=0.0141*
MannW.APAR.Ofra <- wilcox.test(APAR.Sep$Ofra.SH, APAR.Sep$Ofra.DP, paired=FALSE, exact=FALSE)
MannW.APAR.Ofra #APAR is NOT significantly different among sites for Ofra, p=0.2301

#Pigment specific absorption coefficient, a*Chla
a.Sep <- read.xlsx("AllPhysUSVI_T1.xlsx", sheetName="a.")
MannW.a.Sh <- wilcox.test(a.Sep$Ofav.SH, a.Sep$Ofra.SH, paired=FALSE, exact=FALSE)
MannW.a.Sh #a* is significantly different among species in the shallow site, p=0.006167**
MannW.a.Dp <- wilcox.test(a.Sep$Ofav.DP, a.Sep$Ofra.DP, paired=FALSE, exact=FALSE)
MannW.a.Dp #a* is significantly different among species in the deep site, p=0.005817**
MannW.a.Ofav <- wilcox.test(a.Sep$Ofav.SH, a.Sep$Ofav.DP, paired=FALSE, exact=FALSE)
MannW.a.Ofav #a* is NOT significantly different among sites for Ofav, p=0.3722
MannW.a.Ofra <- wilcox.test(a.Sep$Ofra.SH, a.Sep$Ofra.DP, paired=FALSE, exact=FALSE)
MannW.a.Ofra #a* is NOT significantly different among sites for Ofra, p=0.5994


##REGRESSION MODELS TO EXPLORE RELATIONSHIPS OF CHL FLUORESCENCE WITH DEPTH##
#Maximum quantum yield of PSII, Fv/Fm
InteractionFvFm = lm(Yield ~ Depth*Spp, data=FvFmGrad); summary(InteractionFvFm)
#The slopes were not significantly different t(50)=-1.07, p=0.2903
#The intercepts were significantly different t(50)=72.72, p<0.001***
SlopesFvFm = lstrends(InteractionFvFm, "Spp", var="Depth"); SlopesFvFm

#Effective quantum yield of PSII at noon, DF/Fm'
Interaction_dFFm = lm(Yield ~ Depth*Spp, data=dFFmGrad); summary(Interaction_dFFm)
#The slopes were significantly different t(48)=3.68, p=0.000586***
#The intercepts were significantly different t(48)=15.49, p<0.001***
Slopes_dFFm = lstrends(Interaction_dFFm, "Spp", var="Depth"); Slopes_dFFm

#Maximum excitation pressure over PSII at noon, Qm
Interaction_Qm = lm(Yield ~ Depth*Spp, data=QmGrad); summary(Interaction_Qm)
#The slopes were significantly different t(48)=-3.85, p=0.000345***
#The intercepts were significantly different t(48)=10.01, p<0.001***
Slopes_Qm = lstrends(Interaction_Qm, "Spp", var="Depth"); Slopes_Qm


##REGRESSION MODELS TO EXPLORE RELATIONSHIPS OF CHL CONTENT AND OPTICAL PROPERTIES##
#Relationship between Chla content per area and pigment specific absorption coefficient (a*Chla)
a.ALL <- read.xlsx("AllPhysUSVI_T1.xlsx", sheetName="a.ALL") #Analysis for all samples
Interaction.a.ALL = lm(log(a) ~ Chla*Spp, data=a.ALL); summary(Interaction.a.ALL)
#The slopes were not significantly different t(41)=-1.52, p=0.1367
#The intercepts were significantly different t(41)=-17.16, p<0.001***
coefficients(Interaction.a.ALL)
Slopes.a.ALL = lstrends(Interaction.a.ALL, "Spp", var="Chla"); Slopes.a.ALL

a.SH <- read.xlsx("AllPhysUSVI_T1.xlsx", sheetName="a.SH") #Analysis for the shallow site
Interaction.a.SH = lm(log(a) ~ Chla*Spp, data=a.SH); summary(Interaction.a.SH)
#The slopes were not significantly different t(19)=-1.22, p=0.23611
#The intercepts were significantly different t(19)=-9.44, p<0.001***
Slopes.a.SH = lstrends(Interaction.a.SH, "Spp", var="Chla"); Slopes.a.SH

a.DP <- read.xlsx("AllPhysUSVI_T1.xlsx", sheetName="a.DP") #Analysis for the deep site
Interaction.a.DP = lm(log(a) ~ Chla*Spp, data=a.DP); summary(Interaction.a.DP)
#The slopes were not significantly different t(18)=-0.41, p=0.688
#The intercepts were significantly different t(18)=-12.07, p<0.001***
Slopes.a.DP = lstrends(Interaction.a.DP, "Spp", var="Chla"); Slopes.a.DP

#Relationship between Chla content per symbiont and pigment specific absorption coefficient (a*Chla)
a.Sym.ALL <- read.xlsx("AllPhysUSVI_T1.xlsx", sheetName="a.Sym.ALL") #Analysis for all samples
Interaction.a.Sym.ALL = lm(log(a) ~ Ci*Spp, data=a.Sym.ALL); summary(Interaction.a.Sym.ALL)
#The slopes were not significantly different t(38)=-0.404, p=0.6887
#The intercepts were significantly different t(38)=-12.722, p<0.001***
Slopes.a.Sym.ALL = lstrends(Interaction.a.Sym.ALL, "Spp", var="Ci"); Slopes.a.Sym.ALL

a.Sym.SH <- read.xlsx("AllPhysUSVI_T1.xlsx", sheetName="a.Sym.SH") #Analysis for the shallow site
Interaction.a.Sym.SH = lm(log(a) ~ Ci*Spp, data=a.Sym.SH); summary(Interaction.a.Sym.SH)
#The slopes were not significantly different t(18)=0.78, p=0.4431
#The intercepts were significantly different t(18)=-7.153, p<0.001***
Slopes.a.Sym.SH = lstrends(Interaction.a.Sym.SH, "Spp", var="Ci"); Slopes.a.Sym.SH

a.Sym.DP <- read.xlsx("AllPhysUSVI_T1.xlsx", sheetName="a.Sym.DP")
Interaction.a.Sym.DP = lm(log(a) ~ Ci*Spp, data=a.Sym.DP); summary(Interaction.a.Sym.DP)
#The slopes were not significantly different t(16)=-1.127, p=0.2763
#The intercepts were significantly different t(16)=-5.496, p<0.001***
Slopes.a.Sym.SH = lstrends(Interaction.a.Sym.SH, "Spp", var="Ci"); Slopes.a.Sym.SH
