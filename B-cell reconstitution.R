############################
#loading dataset
############################
df <- read.csv(url())

head(df)
summary(df)
dim(df)

#number of repeated measurements per patient
Length <- function (x) sum(x)
ns <- with(df, tapply(bcell, ID, length))
ns <- ns[!is.na(ns)]
c(summary(ns), sd = sd(ns))
boxplot(ns, col = "lightgrey", ylab = "Number Repeated Measurements per Patient")

#plot 16 patients with more than 5 observations
set.seed(321)
ids <- sample(names(ns)[ns >= 5],16) 
xyplot(bcell ~ time| ID,
       panel = function (x, y, ...) {
         panel.xyplot(x, y, type = "p", col = 1, ...)
         panel.loess(x, y, col = 2, lwd = 2, span = 0.8)
       }, data = df, subset = ID %in% ids, layout = c(4, 4), as.table = TRUE, 
       xlab = "Time (days)", ylab = "B cell per ?l")

#mean follow up time B cell
futperid <- with(df, tapply(time, ID, max))

fut <- mean(futperid)

#mean follow up time MB cell
futperid_2 <- with(dfm, tapply(time, ID, max))

fut_2 <- mean(futperid_2[!is.na(futperid_2)])

#plot evolutions
ggplot(df, aes(x = time, y= bcell)) + 
  geom_line(size = 1, alpha = 0.5, aes(shape = ID)) +
  geom_smooth(method = "lm",formula= y ~ ns(x, 3), size = 2, colour="red")  +
  scale_x_continuous(name = "Follow-up time (years)") +
  scale_y_continuous(name = "Total B (cells/µl)",limits = c(0,5000)) +
  theme_bw()

ggplot(df[df$memory<245,], aes(x = time, y= memory)) + 
  geom_line(size = 1, aes(shape = ID)) +
  geom_smooth(method = "lm",formula= y ~ ns(x, 3), size = 2, colour="red")  +
  scale_x_continuous(name = "Follow-up time (years)") +
  scale_y_continuous(name = "Memory B (cells/µl)",limits = c(0,250)) +
  theme_bw()

#plots of B and memory B evolution
library(ggpubr)
ggarrange(p1,p2, pp1,pp2, labels = c("A","B", "C", "D"))
