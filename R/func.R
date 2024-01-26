library(data.table)
library(lme4)
library(ggplot2)
set.seed(1354)
##Simulate population of crosses
addVar <- 5.5
mu <- 33
domVar <- 3
residVar <- 2
numPar <- 20
numCross <- 200
a <- rnorm(n = numPar, 0, sd = sqrt(addVar)) # additive effects
d <- rnorm(n = numCross, 0, sd = sqrt(domVar)) # dominance effects (only in crosses)
#Assuming same residual variance for progeny and parents. Just separating for simplicity
e1 <- rnorm(n = numPar, 0, sd = sqrt(residVar))
e2 <- rnorm(n = numCross, 0, sd = sqrt(residVar))
#make crosses
cr <- as.data.table(expand.grid(1:20, 1:20))
# Random selection of crosses. Non-selfing
crList <- cr[Var1 != Var2][sample(.N, size = numCross)]
setnames(crList, c("Par1", "Par2"))
crList[, y_par1 := 2 * a[Par1] + e1[Par1] + mu
       ][, y_par2 := 2 * a[Par2] + e1[Par2] + mu
        ][, y_progeny := a[Par1] + a[Par2] + d + e2 + mu
          ][, mpv := (y_par1 + y_par2) / 2]
# Regress progney on mid-parent value
ggplot(crList, aes(x = mpv, y = y_progeny)) +
  geom_point() + theme_minimal() + xlab("Mid-Parent Value") +
  ylab("Progeny Value")
lm(y_progeny ~ mpv, data = crList) # Estimated heritability
2 * addVar / (2 * addVar + domVar + residVar) #True heritability


#Instead of these being real progeny values, let's assume they are plot means
#The residual variance now can be interpreted as a within-plot error
# Let's expand the dataset to add a replicate
crList[, cross_id := paste0(letters[Par1], letters[Par2])
       ][, c("Par1", "Par2") := lapply(.(Par1, Par2), as.factor)]
dat <- rbind(crList, crList)
dat[, block := as.factor(rep(1:2, each = 200))
    ][, y_progeny := a[Par1] + a[Par2] + d + mu + rnorm(n = .N, 0, sd = sqrt(residVar))
      ][block == 1, y_progeny := y_progeny + 1.3] #block effect

##Begin with simple model
mod1 <- lmer(y_progeny ~ block + (1 | cross_id), data = dat)
tmp <- as.data.frame(VarCorr(mod1))$vcov
tmp[1] / sum(tmp) #Estimated broad-sense heritability
(2 * addVar + domVar) / (2 * addVar + domVar + residVar) #True broad-sense heritability

##We can also leverage parental information to get narrow-sense
mod2 <- lmer(y_progeny ~ block + (1 | cross_id) + (1 | Par1) + (1 | Par2), data = dat)
summary(mod2)
##Before, we assumed additive effects with the same variance for each parent.
#This is harder without specifying model matrices explicitly in lmer
#Assuming that var_add = 2*var_gca (different in outbred tetras)
## See (Endelman et al, 2018)
tmp2 <- as.data.frame(VarCorr(mod2))$vcov
sum(tmp2[2:3]) / (sum(tmp2[-4]) + tmp2[4] / 2) #Estimated narrow-sense heritability
(2 * addVar) / (2 * addVar + domVar + residVar) #True narrow-sense heritability
