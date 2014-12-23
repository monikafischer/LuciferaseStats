#command-shift-o collapses everything!
# t.tests and anovas on samples within a dataframe using rbind #####
LucPrelimDat <- read.table(file.choose(), header=TRUE)

WT <- LucPrelimDat[LucPrelimDat$strain=="WT9",]
ham11 <- LucPrelimDat[LucPrelimDat$strain=="ham11_3.12.10",]
WTvHAM11 <- rbind(WT, ham11)
t.test(cps~strain, data=WTvHAM11)
#p=1.606e-5

lucdat <- read.table(file.choose(), header=TRUE)

h314_ham11 <- lucdat[lucdat$strain=="h314_ham11",]
h314_2489 <- lucdat[lucdat$strain=="h313_2489",]
h312_ham11 <- lucdat[lucdat$strain=="h312_ham11",]
h312_2489 <- lucdat[lucdat$strain=="h312_2489",]
WT9_ham11 <- lucdat[lucdat$strain=="WT9_ham11",]
WT9_2489 <- lucdat[lucdat$strain=="WT9_2489",]
WT7_ham11 <- lucdat[lucdat$strain=="WT7_ham11",]
WT7_2489 <- lucdat[lucdat$strain=="WT7_2489",]

h314 <- rbind(h314_ham11, h314_2489)
t.test(cps~strain, data=h314)

##ANOVAs
mix <- rbind(WT7_ham11, WT9_ham11, h314_2489, h312_2489)
anova <- aov(cps~strain, data=mix)
summary(anova)
#p=0.0473
TukeyHSD(anova)

lumdat <- read.table(file.choose(), header=TRUE)
head(lumdat)
anova <- aov(luc~strain, data=lumdat)
summary(anova)
TukeyHSD(anova)

#Reshape data directly from the luminometer: ####
#each collumn is a strain, each row is a rep ####
library("reshape2", lib.loc="/Library/Frameworks/R.framework/Versions/3.0/Resources/library")
dat <- read.table(file.choose(), header=TRUE)
head(dat)
datm <- melt(dat, id.vars="vars")
head(datm)
#collumn 1 is called "vars" and was added in Excel. It's just the numbers 1-8.
anova <- aov(value~variable, data=datm)
summary(anova)
tukey <- TukeyHSD(anova, "variable")
result <- data.frame(tukey$variable)
result["p.adj"]
#prints only the important stuff from the Tukey test!

nrc1.3 <- datm[datm$variable=="nrc1.3",]
bem1.4 <- datm[datm$variable=="bem1.4",]
nrc1bem1 <- rbind(nrc1.3, bem1.4)
t.test(value~variable, data=nrc1bem1)
t.test(value~variable, data=datm[which(datm$variable %in% c("nrc1.3", "bem1.4")), ])

#Alex's two cents of awesome:
#Datatables!
#with data.table you can avoid creating a million objects with things like rbind
datm.dt<-data.table(datm)
t.test(value~variable, data=datm.dt[variable %in% c("nrc1.3", "bem1.4"), ])
datm.dt
datm.dt[ ,sd(value)/sqrt(.N), by=variable][order(V1, decreasing=TRUE), ]
datm.dt[ ,.N, by=variable]
datm.mean.se<-datm.dt[ ,list(mean(value), 
                             mean(value)+sd(value/sqrt(.N)), 
                             mean(value)-sd(value/sqrt(.N))), by=variable]

theme_opts <- list(theme(panel.background = element_rect(fill="white"),
                          panel.grid.minor = element_blank(),
                          panel.grid.major = element_blank(),
                          axis.title.x = element_blank()))
                  
ggplot(datm.mean.se)+
  geom_bar(aes(variable, V1), fill="orange", color="red", stat="identity")+
  geom_errorbar(aes(x=variable, ymax=V2, ymin=V3), color="purple", width=0.3)+
  theme_opts+
  ylab("Pprm-1::Luciferase Output (photons/sec)")

library(plotrix)
gap.barplot(datm.mean.se$V1, gap=c(400,2700))
#broken y-axis ...this can't be done in ggplot2

#each row is a strain, each collumn is a rep #######
dat <- read.table(file.choose(), header=TRUE)
head(dat)
#rotate table so that each collumn is a strain and each row is a rep:
datcast <- dcast(melt(dat), variable~strain)
head(datcast)
#after using dcast the first row is the vars, remove it:
dat2 <- datcast[-1,]
head(dat2)
#dcast also renamed the first collumn to "variable", which gets confusing later with melt,
#so I changed the name of the first collumn back to "vars":
colnames(dat2)[1] <- "vars"
datm <- melt(dat2, id.vars="vars")
head(datm)
anova <- aov(value~variable, data=datm)
summary(anova)
TukeyHSD(anova)
tukey <- TukeyHSD(anova)
head(tukey)
result <- data.frame(tukey$variable)
head(result)
result["p.adj"]
#something about printing only the p values changes the way the
#p-values are formated, for example 0.0000146 --> 1.457372e-05

# t.tests on a dataframe made by mannually entering in values with c() and rbind ##########

#This code gives a different answer and I'm not sure why....
#WT <- c(4520, 4600, 2720, 3760, 3320, 5320, 4000, 3240, 9600, 12160, 11640)
#ham11 <- c(680, 600, 560, 440, 560, 640, 840, 440, 1880, 2200, 1800)
#t.test(WT, ham11)

#This seems like the best way to do it...
#data.frame recognizes that it's what I want it to be
WT <- c(4520, 4600, 2720, 3760, 3320, 5320, 4000, 3240, 9600, 12160, 11640)
ham11 <- c(680, 600, 560, 440, 560, 640, 840, 440, 1880, 2200, 1800)
ham1310 <- c(1440, 1560, 1240, 1160, 1120, 1520, 1280, 1000, 2240, 2200, 2080)
ham139 <- c(1320, 1080, 920, 1080, 1000, 960, 1080, 920, 920, 1440, 1080)
lum <- c(ham11, ham1310)
strain <- c("ham11", "ham1310")
strain <- rep(strain, c(11,11))
dataframe <- data.frame(lum, strain)
t.test(lum~strain, data=dataframe)

#This works, but it seems messier
#matrix seems very barebones, it's just a matrix
strain2 <- c("WT", "WT", "WT", "WT", "WT", "WT", "WT", "WT", "WT", "WT", "WT", "ham11", "ham11", "ham11", "ham11", "ham11", "ham11", "ham11", "ham11", "ham11", "ham11", "ham11")
lum2 <- c(4520, 4600, 2720, 3760, 3320, 5320, 4000, 3240, 9600, 12160, 11640, 680, 600, 560, 440, 560, 640, 840, 440, 1880, 2200, 1800)
dat <- cbind (strain2, lum2)
datmatrix <- matrix(data = dat, nrow = 22, ncol = 2)
t.test(lum2~strain2, data=datmatrix)
