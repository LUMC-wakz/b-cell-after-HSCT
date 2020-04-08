#########################################
#plots at 1 year
#########################################
#packages
library(ggplot2)
library(dplyr)
library(ggpubr)

########################################
#data at 1 year
########################################
df12m <- df %>% filter(between(time, 0.8,1.2)) %>% distinct(ID, .keep_all = T)

######################################
#B cells
######################################

#plot differences in age donor 
p1 <- ggplot(df12m, aes(x = donorage, y = bcell, fill = donorage)) + geom_boxplot() + 
  scale_x_discrete(name = "Donor age") +
  scale_y_continuous(name = "Cells/µl", limits = c(0,5000)) +
  theme_bw() + stat_compare_means(label.x = 1.2) + theme(legend.position = 'none',
                                            title = element_text(size=24,face="bold"),
                                            axis.text=element_text(size=12),
                                            axis.title=element_text(size=14,face="bold"))

#age recipient
p2 <- ggplot(df12m, aes(x = recage, y = bcell, fill = recage)) + geom_boxplot() + 
  scale_x_discrete(name = "Recipient age") +
  scale_y_continuous(name = "Cells/µl", limits = c(0,5000)) +
  theme_bw() + stat_compare_means(label.x = 1.2) + theme(legend.position = 'none',
                                                                       title = element_text(size=24,face="bold"),
                                                                       axis.text=element_text(size=12),
                                                                       axis.title=element_text(size=14,face="bold"))

#plot differences in source 
p3 <- ggplot(df12m, aes(x = Source, y = bcell, , fill = Source)) + geom_boxplot() + 
  scale_x_discrete(name = "Source", labels = c("BM","CB", "PBSC")) +
  scale_y_continuous(name = "Cells/µl", limits = c(0,5000)) +
  theme_bw() + stat_compare_means() + theme(legend.position = 'none',
                                           title = element_text(size=24,face="bold"),
                                           axis.text=element_text(size=12),
                                           axis.title=element_text(size=14,face="bold"))

#differences in dontype
p4 <- ggplot(df12m, aes(x = Dontype, y = bcell, fill = Dontype)) + geom_boxplot() + 
  scale_x_discrete(name = "Donor type", labels = c("Haplo","IRD", "MUD")) +
  scale_y_continuous(name = "Cells/µl", limits = c(0,5000)) +
  theme_bw() + stat_compare_means() + theme(legend.position = 'none',
                                                                      title = element_text(size=24,face="bold"),
                                                                      axis.text=element_text(size=12),
                                                                      axis.title=element_text(size=14,face="bold"))

#conditioning
p5 <- ggplot(df12m, aes(x = con, y = bcell, fill = con)) + geom_boxplot() + 
  scale_x_discrete(name = "Conditioning") +
  scale_y_continuous(name = "Cells/µl", limits = c(0,5000)) +
  theme_bw() + stat_compare_means() + theme(legend.position = 'none',  title = element_text(size=24,face="bold"),
                                                                       axis.text=element_text(size=12),
                                                                       axis.title=element_text(size=14,face="bold"))


######################################
#memory B cells
######################################

#plot differences in age donor 
p6 <- ggplot(df12m, aes(x = donorage, y = memory, fill = donorage)) + geom_boxplot() + 
  scale_x_discrete(name = "Donor age") +
  scale_y_continuous(name = "Cells/µl", limits = c(0,250)) +
  theme_bw() + stat_compare_means(label.x = 1.2) + theme(legend.position = 'none',
                                            title = element_text(size=24,face="bold"),
                                            axis.text=element_text(size=12),
                                            axis.title=element_text(size=14,face="bold"))

#age recipient
p7 <- ggplot(df12m, aes(x = recage, y = memory, fill = recage)) + geom_boxplot() + 
  scale_x_discrete(name = "Recipient age") +
  scale_y_continuous(name = "Cells/µl", limits = c(0,250)) +
  theme_bw() + stat_compare_means(label.x = 1.2) + theme(legend.position = 'none',
                                            title = element_text(size=24,face="bold"),
                                            axis.text=element_text(size=12),
                                            axis.title=element_text(size=14,face="bold"))

#plot differences in source 
p8 <- ggplot(df12m, aes(x = Source, y = memory, fill = Source)) + geom_boxplot() + 
  scale_x_discrete(name = "Source", labels = c("BM","CB", "PBSC")) +
  scale_y_continuous(name = "Cells/µl", limits = c(0,250)) +
  theme_bw() + stat_compare_means() + theme(legend.position = 'none',
                                            title = element_text(size=24,face="bold"),
                                            axis.text=element_text(size=12),
                                            axis.title=element_text(size=14,face="bold"))

#differences in dontype
p9 <- ggplot(df12m, aes(x = Dontype, y = memory, fill = Dontype)) + geom_boxplot() + 
  scale_x_discrete(name = "Donor type", labels = c("Haplo","IRD", "MUD")) +
  scale_y_continuous(name = "Cells/µl", limits = c(0,250)) +
  theme_bw() + stat_compare_means() + theme(legend.position = 'none',
                                            title = element_text(size=24,face="bold"),
                                            axis.text=element_text(size=12),
                                            axis.title=element_text(size=14,face="bold"))

#conditioning
p10 <- ggplot(df12m, aes(x = con, y = memory, fill = con)) + geom_boxplot() + 
  scale_x_discrete(name = "Conditioning") +
  scale_y_continuous(name = "Cells/µl", limits = c(0,250)) +
  theme_bw() + stat_compare_means() + theme(legend.position = 'none',
                                                                       title = element_text(size=24,face="bold"),
                                                                       axis.text=element_text(size=12),
                                                                       axis.title=element_text(size=14,face="bold"))

#arrange them in one plot
Fig_2AJ <- ggarrange(p1,p6,p2,p7,p3,p8,p4,p9,p5,p10, labels = c("A","B","C","D","E","F","G","H","I","J"), ncol = 2, nrow = 5)
Fig_2AJ