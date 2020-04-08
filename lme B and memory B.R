##################################
#model fitting 0.5-2 years
##################################
#packages
library("lattice")
library("nlme")
library("splines")
library("ggplot2")
library("Hmisc")
library("dplyr")

#load data
df <- read.csv(url())

#plot evolutions
p1 <- ggplot(df, aes(x = time, y= bcell)) + 
  geom_line(size = 1, alpha = 0.5, aes(shape = ID)) +
  geom_smooth(method = "lm",formula= y ~ ns(x, 3), size = 2, colour="red")  +
  scale_x_continuous(name = "Follow-up time (years)") +
  scale_y_continuous(name = "total B (cells/µl)",limits = c(0,5000)) +
  theme_bw()

p2 <- ggplot(df[df$memory<245,], aes(x = time, y= memory)) + 
  geom_line(size = 1, alpha = 0.5, aes(shape = ID)) +
  geom_smooth(method = "lm",formula= y ~ ns(x, 3), size = 2, colour="red")  +
  scale_x_continuous(name = "Follow-up time (years)") +
  scale_y_continuous(name = "memory B (cells/µl)",limits = c(0,250)) +
  theme_bw()

Fig_3AB <- ggarrange(p1,p2, labels = c("A","B"), ncol = 2, nrow = 1)
Fig_3AB

#fitting a lme
fm_1 <- lme(sqrt(bcell) ~ ns(time, 3) + AgeD + Source + Dontype + AgeR + con, data = df, 
            random = ~ 1 | ID, control=lmeControl(opt = "optim"))

summary(fm_1)

fm_2 <- update(fm_1, random = ~ time | ID) #Time random slopes

anova(fm_1, fm_2) #significant difference

fm_3 <- update(fm_1, random = ~ ns(time, 2) | ID) #add 2nd degree spline to random slopes

anova(fm_2, fm_3) #significant difference

fm_4 <- update(fm_1, random = ~ ns(time, 3) | ID) #add 3nd degree spline to random slopes

anova(fm_3, fm_4) #not significant, proceed with fm_3

#fit under ML instead of REML to do likelihood ratio test
fm_3ML <- update(fm_3, method = "ML")

#now we test the nonlinear effects
fm_3ML_Lin <- update(fm_3ML, fixed = . ~ time + AgeD  + Source + Dontype + AgeR)

anova(fm_3ML, fm_3ML_Lin) #sign

#proceed with fm_3ML

summary(fm_3ML)

#fixed effects with 95%CI
intervals(fm_3ML)

#effect plots for the model

effectPlotData <- function (object, newdata, orig_data) {
  form <- formula(object)
  namesVars <- all.vars(form)
  betas <- if (!inherits(object, "lme")) coef(object) else fixef(object)
  V <- if (inherits(object, "geeglm")) object$geese$vbeta else vcov(object)
  orig_data <- orig_data[complete.cases(orig_data[namesVars]), ]
  Terms <- delete.response(terms(form))
  mfX <- model.frame(Terms, data = orig_data)
  Terms_new <- attr(mfX, "terms")
  mfX_new <- model.frame(Terms_new, newdata, xlev = .getXlevels(Terms, mfX))
  X <- model.matrix(Terms_new, mfX_new)
  pred <- c(X %*% betas)
  ses <- sqrt(diag(X %*% V %*% t(X)))
  newdata$pred <- pred
  newdata$low <- pred - 1.96 * ses
  newdata$upp <- pred + 1.96 * ses
  newdata
}

#effect plot
newDF_1 <- with(df2n, expand.grid(time = seq(0.5,2, length.out = 25), 
                                  Source = "Bone Marrow", Dontype = c("Matched unrelated", "Identical related"), 
                                  con = "MAC", AgeR = 8))

newDF_1$AgeD[1:25] = 30
newDF_1$AgeD[26:50] = 5

predictions <- effectPlotData(fm_3ML, newDF_1, df2n)

#ggplot
ep1 <- ggplot(data = predictions, aes(x = time, y = pred ^ 2, col = factor(AgeD))) + 
  geom_line(size = 2) +
  geom_ribbon(data = predictions, aes(ymin = low^ 2, ymax = upp^ 2), linetype = 0, alpha = 0.3) + 
  scale_x_continuous(name = "Follow-up time (years)") +
  scale_y_continuous(name = "Total B (cells/µl)", limits = c(0,1500)) +
  theme_bw() + theme(legend.position = "none") +
  scale_color_manual(values = c("#D95F02", "#7570B3"))


# residuals and qqplots
plot(fm_3ML, resid(., type = "p") ~ fitted(.)| cut(AgeD, breaks = 4), 
     type = c("p", "smooth"), lwd = 3)
plot(fm_3ML, resid(., type = "p") ~ fitted(.)| cut(AgeR, breaks = 4), 
     type = c("p", "smooth"), lwd = 3)
plot(fm_3ML, resid(., type = "p") ~ fitted(.)| con, 
     type = c("p", "smooth"), lwd = 3)
plot(fm_3ML, resid(., type = "p", level = 0) ~ fitted(., level = 0)| cut(AgeD, breaks = 4, labels = c("0-14","14-28", "28-43", "43-60")), 
     type = c("p", "smooth"), lwd = 3)
plot(fm_3ML, resid(., type = "p", level = 0) ~ fitted(., level = 0)| cut(AgeR, breaks = 4, labels = c("0-14","14-28", "28-43", "43-60")), 
     type = c("p", "smooth"), lwd = 3)
plot(fm_3ML, resid(., type = "p", level = 0) ~ fitted(., level = 0)| con, 
     type = c("p", "smooth"), lwd = 3)
qqnorm(fm_3ML, ~ resid(., type = "p"), main = "QQ-plot of the standardized conditional residuals" )
qqnorm(fm_3ML, ~ resid(., type = "p", level = 0), main = "QQ-plot of the standardized marginal residuals " )

#correlation matrices
getVarCov(fm_3ML_no_diagn_interaction, individuals = 175, type = "marginal")
cov2cor(getVarCov(fm_3ML_no_diagn_interaction, individuals = 175, type = "marginal")[[1]])

#########################################################
#memory B cell
#########################################################
#remove missing values
dfm <- df[!is.na(df$memory),]

#fitting a lme
fm_1b <- lme(log(memory) ~ ns(time,3) + AgeD + Source + Dontype + AgeR + con, data = dfm, 
            random = ~ 1 | ID, control=lmeControl(opt = "optim"))

summary(fm_1b)

fm_2b <- update(fm_1b, random = ~ time | ID) #Time random slopes

anova(fm_1b, fm_2b) #significant difference

fm_3b <- update(fm_1b, random = ~ ns(time, 2) | ID) #add 2nd degree spline to random slopes

anova(fm_2b, fm_3b) #significant difference

fm_4b <- update(fm_1b, random = ~ ns(time, 3) | ID) #add 2nd degree spline to random slopes

anova(fm_3b, fm_4b) #significant difference, continue with fm_4b

#fit under ML instead of REML to do likelihood ratio test
fm_4MLb <- update(fm_4b, method = "ML")

#now we test the nonlinear effects
fm_4ML_Linb <- update(fm_4MLb, fixed = . ~ time + Source + AgeR + AgeD + Dontype)

anova(fm_4MLb, fm_4ML_Linb) #sign

#proceed with fm_4MLb

summary(fm_4MLb)

#fixed effects with 95%CI
intervals(fm_4MLb)

#marginal and subject specific prediction from model

predictions2 <- effectPlotData(fm_4MLb, newDF_1, df2n)

#ggplot
ep2 <- ggplot(data = predictions2, aes(x = time, y = exp(pred), col = factor(AgeD))) + 
  geom_line(size = 2) +
  geom_ribbon(data = predictions2, aes(ymin = exp(low), ymax = exp(upp)), linetype = 0, alpha = 0.3) + 
  scale_x_continuous(name = "Follow-up time (years)") +
  scale_y_continuous(name = "Memory B (cells/Âµl)", limits = c(0,150)) +
  theme_bw() + theme(legend.position = "none") +
  scale_color_manual(values = c("#D95F02", "#7570B3"))

# residuals and qqplots
plot(fm_4MLb, resid(., type = "p") ~ fitted(.)| cut(AgeD, breaks = 4), 
     type = c("p", "smooth"), lwd = 3)
plot(fm_4MLb, resid(., type = "p") ~ fitted(.)| cut(AgeR, breaks = 4), 
     type = c("p", "smooth"), lwd = 3)
plot(fm_4MLb, resid(., type = "p") ~ fitted(.)| con, 
     type = c("p", "smooth"), lwd = 3)
plot(fm_4MLb, resid(., type = "p", level = 0) ~ fitted(., level = 0)| cut(AgeD, breaks = 4, labels = c("0-14","14-28", "28-43", "43-60")), 
     type = c("p", "smooth"), lwd = 3)
plot(fm_4MLb, resid(., type = "p", level = 0) ~ fitted(., level = 0)| cut(AgeR, breaks = 4, labels = c("0-14","14-28", "28-43", "43-60")), 
     type = c("p", "smooth"), lwd = 3)
plot(fm_4MLb, resid(., type = "p", level = 0) ~ fitted(., level = 0)| con, 
     type = c("p", "smooth"), lwd = 3)
qqnorm(fm_4MLb, ~ resid(., type = "p"), main = "QQ-plot of the standardized conditional residuals" )
qqnorm(fm_4MLb, ~ resid(., type = "p", level = 0), main = "QQ-plot of the standardized marginal residuals " )

#arrange multiple plots in one plot
library(ggpubr)
#evolutions
Fig_4AB <- ggarrange(ep1,ep2, labels = c("A","B"), ncol = 2, nrow = 1)
Fig_4AB