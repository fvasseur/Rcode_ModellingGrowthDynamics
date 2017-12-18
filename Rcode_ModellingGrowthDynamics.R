
#################################################################################################
#
# Model growth dynamics (M, GR, and RGR) with predicted and measured rosette dry mass
#
#################################################################################################



# Import phenotypic data: estimated dry mass during ontogeny and measured final dry mass (combined in one file)
#-----------------------------------------------------------------------------------------

totflo <- read.csv("DM_ontogeny_AllIndividuals.csv", header=T, sep=",")
totflo$idPot <- as.factor(as.character(totflo$idPot))
totflo$daA02 <- as.factor(as.character(totflo$daA0))
totflo $ros_DM <- 1000 * totflo $ros_DM # to convert in mg
totflo $agDM <- 1000 * totflo $agDM  # to convert in mg
totflo3 <- totflo

tab_1 <- totflo[totflo$daA0==totflo$Growth_duration,]
tab_1 <- tab_1[,c(1:11)]
tab_1 <- droplevels(tab_1)




# Sigmoid growth curve M = f(t)
#................................
pdf("Fitted_Growth_Curves.pdf")
for(j in levels(totflo3$idPot))  
{   
  par(mar=c(.5,.5,.5,.5), oma=c(5,5,.5,.5))
  plot(totflo3 $daA0, totflo3 $ros_DM, type="n", xlab="", ylab="", 
       xlim=c(0,1.1*max(totflo3[totflo3 $idPot %in% c(j),"daA0"], na.rm=T)), 
       ylim=c(0,1.1*max(totflo3[totflo3 $idPot %in% c(j),"ros_DM"], na.rm=T))) #
  mtext("Rosette DM (ng)", side = 2, line = 2.5, outer = T, at = 0.5, cex=1.3)
  mtext("Age (d)", side = 1, line = 2.5, outer = T, at = 0.5, cex=1.3)
  legend("topleft", legend=c(paste("pot", j)), cex=1.5, bty="n", lwd=1, pch=16, col=c("black"))
  
  DMmax <- max(totflo3[totflo3 $idPot %in% c(j),"ros_DM"]) 
  color <- "black"  
  S.mod <- nls(ros_DM ~ DMmax / (1 + exp(-(daA0 - X0) / B)), 
               data = totflo3[ totflo3 $idPot %in% c(j) ,], 
               start = list(B=1, X0=mean(totflo3[totflo3 $idPot %in% c(j), "daA0"])), 
               nls.control(maxiter=1000, minFactor = 1/1000))
  B.sigm <- coef(S.mod)[1]
  X0.sigm <- coef(S.mod)[2]
  tab_1[tab_1$idPot %in% c(j),"ModDM_coefB"] <- coef(S.mod)[1]
  tab_1[tab_1$idPot %in% c(j),"ModDM_coefB_ciL"] <- as.numeric(confint(S.mod)[1,1])
  tab_1[tab_1$idPot %in% c(j),"ModDM_coefB_ciU"] <- as.numeric(confint(S.mod)[1,2])
  tab_1[tab_1$idPot %in% c(j),"ModDM_coefX0"] <- coef(S.mod)[2] 
  tab_1[tab_1$idPot %in% c(j),"ModDM_coefX0_ciL"] <- as.numeric(confint(S.mod)[2,1])
  tab_1[tab_1$idPot %in% c(j),"ModDM_coefX0_ciU"] <- as.numeric(confint(S.mod)[2,2])
  tab_1[tab_1$idPot %in% c(j),"ModDM_Pvalue"] <- S.mod$convInfo$finTol
  tab_1[tab_1$idPot %in% c(j),"ModDM_AIC"] <- AIC(S.mod)
  curve.S <- function(age, A, B, X0) { A / (1 + exp(-(age - X0) / B)) }
  curve(curve.S(age = x, A = DMmax, B = B.sigm, X0 = X0.sigm), 0, 
        max(totflo3[totflo3 $idPot %in% c(j),"daA0"]), lwd = 1, add=T, col=color)
  points(totflo3[totflo3 $idPot==j, "daA0"], totflo3[totflo3 $idPot==j, "ros_DM"], 
         type="p", pch=16, cex=0.5, col=color)
  print(j)
}
dev.off()

dat_GrowthParam_GT01 <- tab_1

write.table(dat_GrowthParam_GT01,"dat_GrowthParam_GT01.csv", sep=",", dec=".")




# -----------------------------------------------------
# Calculation of M(t), GR(t) and RGR(t) for a pot j
# -----------------------------------------------------


# Define pot (here pot ID = 173)
j = 173

# Model sigmoid growth curve and extract parameters
DMmax <- max(totflo3[totflo3 $idPot %in% c(j),"ros_DM"]) 
S.mod <- nls(ros_DM ~ DMmax / (1 + exp(-(daA0 - X0) / B)), 
             data = totflo3[ totflo3 $idPot %in% c(j) ,], 
             start = list(B=1, X0=mean(totflo3[totflo3 $idPot %in% c(j), "daA0"])), 
             nls.control(maxiter=1000, minFactor = 1/100000000000))
B.sigm <- coef(S.mod)[1]
B.sigm_ciL  <- as.numeric(confint(S.mod)[1,1])
B.sigm_ciU <- as.numeric(confint(S.mod)[1,2])
X0.sigm <- coef(S.mod)[2]
X0.sigm_ciL  <- as.numeric(confint(S.mod)[2,1])
X0.sigm_ciU <- as.numeric(confint(S.mod)[2,2])


# time sequence for life cyle duration of pot j
seq1 <- seq(0,max(totflo3[ totflo3 $idPot %in% c(j) ,"daA0"]),1)

# M(t)
tDM <- DMmax / (1 + exp(-(seq1 - X0.sigm) / B.sigm))
tDM_ciL <- DMmax / (1 + exp(-(seq1 - X0.sigm_ciL) / B.sigm_ciL))
tDM_ciU <- DMmax / (1 + exp(-(seq1 - X0.sigm_ciU) / B.sigm_ciU))

# GR(t)
tGR <- (tDM/B.sigm) - ((tDM ^2)/(B.sigm*DMmax))
tGR_ciL <- (tDM_ciL/B.sigm_ciL) - ((tDM_ciL ^2)/(B.sigm_ciL*DMmax))
tGR_ciU <- (tDM_ciU/B.sigm_ciU) - ((tDM_ciU ^2)/(B.sigm_ciU*DMmax))

# RGR(t)
tRGR <- (tGR/tDM)*1000 # to convert RGR in mg d-1 g-1
tRGR_ciL <- (tGR_ciL/tDM_ciL)*1000 # to convert RGR in mg d-1 g-1
tRGR_ciU <- (tGR_ciU/tDM_ciU)*1000 # to convert RGR in mg d-1 g-1
