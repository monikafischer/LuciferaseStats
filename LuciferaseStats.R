#command-option-o collapses everything!
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

#Alex's two cents of awesome: Datatables! ####
#with data.table you can avoid creating a million objects with things like rbind
library("data.table", lib.loc="/Library/Frameworks/R.framework/Versions/3.0/Resources/library")
datm.dt<-data.table(datm)
t.test(value~variable, data=datm.dt[variable %in% c("nrc1.3", "bem1.4"), ])
datm.dt
datm.dt[ ,sd(value)/sqrt(.N), by=variable][order(V1, decreasing=TRUE), ]
datm.dt[ ,.N, by=variable]
datm.mean.se<-datm.dt[ ,list(mean(value), 
                             mean(value)+sd(value/sqrt(.N)), 
                             mean(value)-sd(value/sqrt(.N))), by=variable,]

#datatable basics ####
library("data.table", lib.loc="/Library/Frameworks/R.framework/Versions/3.0/Resources/library")
dat <- read.table(file.choose(), header=TRUE) #read in summary luciferase file
datnovars <- dat[ ,-1] #remove first column, which is "vars"
head(datnovars)
dat.dt <- data.table(datnovars)
dat.dt #automatically shows you the first 5 rows and the last 5 rows! cool!
dat.dt[5] #outputs row #5  ...dat.dt[5,] also works. Commas are optional.
dat.dt[strain=="WT9"] #outputs all rows with "WT9" in the strain column!
dat.dt[date=="6.Feb.14"] #outputs all rows with "6.Feb.14" in the date column
setkey(dat.dt, strain) #sorts entire table based on the strain column!
dat.dt["WT9"] #after setting the key, there's no need to use "strain==" to pull rows from the key column
setkey(dat.dt, strain, date) #re-set the key to sort by two columns instead of one.
dat.dt[list("WT9", "6.Feb.14")] #outputs row with "WT9" in the strain column and "the date"6.Feb.14" in the date column

#
# plotting with ggplot2 ####
theme_opts <- list(theme(panel.background = element_rect(fill="white"),
                          panel.grid.minor = element_blank(),
                          panel.grid.major = element_blank(),
                          axis.title.x = element_blank()))
                          #plot.background = element_blank(),
                          #panel.border = element_blank(),
                          #panel.margin = unit(1, "lines"),
                          #axis.line = element_line(size=0.2, linetype=1, color="black")))
                          #axis.text.x = element_blank(),
                          #axis.text.y = element_blank(),
                          #axis.ticks = element_blank(),
                          #axis.title.x = element_blank(),
                          #axis.title.y = element_blank(),
                          #plot.title = element_text(size=22)))
                          #xlab(expression(paste(" ", Delta, "gene"))) #add Greek Delta to title!

#basic barplot of sample means with error bars
library("ggplot2", lib.loc="/Library/Frameworks/R.framework/Versions/3.0/Resources/library")
ggplot(datm.mean.se)+
  geom_bar(aes(variable, V1), fill="orange", color="red", stat="identity")+
  geom_errorbar(aes(x=variable, ymax=V2, ymin=V3), color="purple", width=0.3)+
  theme_opts+
  ylab("Pprm-1::Luciferase Output (photons/sec)")


#Reshape data directly from the luminometer: ####
#each row is a strain, each collumn is a rep #######
library("reshape2", lib.loc="/Library/Frameworks/R.framework/Versions/3.0/Resources/library")
rm(dat)
dat <- read.table(file.choose(), header=TRUE)
head(dat)
datcast <- dcast(melt(dat), variable~strain) #rotate table so collumns=strains and rows=reps
head(datcast)
dat2 <- datcast[-1,] #after using dcast the first row is the vars, remove it
head(dat2)
#dcast also renamed the first collumn to "variable", which gets confusing later with melt,
#so I changed the name of the first collumn back to "reps":
colnames(dat2)[1] <- "reps"
datm <- melt(dat2, id.vars="reps") #reps is now a column that functions as the vars
head(datm)
anova <- aov(value~variable, data=datm)
summary(anova)
tukey <- TukeyHSD(anova)
result <- data.frame(tukey$variable)
head(result)
result["p.adj"]
#something about printing only the p values changes the way the
#p-values are formated, for example 0.0000146 --> 1.457372e-05 ...but not always...
plot(TukeyHSD(anova, "variable")) #plots TukeyHSD p-values...

###
### ANOVAs on dataframe containing all of my raw luciferase data

totlumdat <- read.table(file.choose(), header=TRUE)
datcast <- dcast(melt(totlumdat[totlumdat$date=="12.Jan.15",]), variable~strain)
dat2 <- datcast[-1,]
colnames(dat2)[1] <- "reps"
datm <- melt(dat2, id.vars="reps") #reps is now a column that functions as the vars
datmo <- na.omit(datm) #remove any row that contains an NA
anova <- aov(value~variable, data=datmo)
summary(anova)
tukey <- TukeyHSD(anova)
result <- data.frame(tukey$variable)
result["p.adj"]
#
datmoNoWT <- datmo[-c(71:80), ]
datmoNoWT
anova <- aov(value~variable, data=datmoNoWT)
summary(anova)
tukey <- TukeyHSD(anova)
result <- data.frame(tukey$variable)
result["p.adj"]
#
datmoWTonly <- datmo[-c(1:66), ]
datmoWTonly
t.test(value~variable, data=datmoWTonly)

#### alex had to play with this -------
alex.dat<-dat.dt
#create ordered factor
alex.dat[ ,group_f:=factor(group, levels=c('low', 'med', 'high'))]
alex.dat
#create plot
alex.p<-ggplot(alex.dat)+
  geom_point(aes(strain, normavg, color=group_f))+
  scale_x_discrete(limits=alex.dat[order(group_f), unique(strain)])
alex.p
#
alex.dat<-dat.dt
#create ordered factor
alex.dat[ ,group_f:=factor(group, levels=c('low', 'med', 'high'))]
alex.dat
#create plot
alex.p<-ggplot(alex.dat)+
  geom_point(aes(strain, normavg, color=group_f))+
  scale_x_discrete(limits=alex.dat[order(group_f), unique(strain)])
alex.p

#saturday morning mess around
head(alex.dat)
alex.dat[ ,lapply(.SD, class)]
summary(alex.dat)
alex.dat[ ,rep12:=NULL]

alex.melt<-melt(alex.dat, id.vars=c("date", "strain", "gene", "isolate", "group", "group_f"), na.rm=TRUE)
alex.melt
ggplot(alex.melt)+
  geom_point(aes(gene, value, color=group_f))+
  scale_x_discrete(limits=alex.melt[order(group_f), unique(gene)])+
  theme(axis.text.x = element_text(angle=90))


#import raw data?
install.packages("XML")
library(XML)
dat.raw<-htmlTreeParse("/Users/monikafischer/Dropbox/GlassLab/Luciferase/Raw_Data_Files/2Dec2013.mht")
dat.raw<-readLines("/Users/monikafischer/Dropbox/GlassLab/Luciferase/Raw_Data_Files/2Dec2013.mht")
dat.raw[grep("</TABLE>", dat.raw)]
#never mind, kinda complicated, check if can export in different format...

## plotting normalized luciferase data ####
normdat <- read.table(file.choose(), header=TRUE)
head(normdat)
normdat.dt <- data.table(normdat)
setkey(normdat.dt, date)
normdat.dt.noaug21 <- normdat.dt[!"21.Aug.14"]
normdat.dt.noaug21may18 <- normdat.dt.noaug21[!"18.May.14"] 
#remove all rows containing 21.Aug.14 and 18.May.14 because these experiments lacked mak2 as a control.
#normdat.dt[ ,rep12:=NULL] #remove column 12 because it's all NA's ...removed in Excel prior to importing
normdat.dt.noaug21may18[ ,group_f:=factor(group, levels=c('low', 'med', 'high'))]
#creates an ordered factor which will ultimately dictate the order inwhich the groups appear on a plot

#boxplot of normalized averages:
plot(normavg~group_f, data=normdat.dt.noaug21may18, 
     ylab="Normalized Average Luciferase Output", xlab="Group")

#barplot of nomalized averages with sd:
library(plyr)
sdmeannormavg <- ddply(normdat.dt.noaug21may18, "group_f", summarise, mnormavg = mean(normavg), sdnormavg = sd(normavg))
ggplot(sdmeannormavg, aes(x=factor(group_f), y=mnormavg))+
  geom_bar(stat = "identity", width = 0.3, fill="white", color="black")+  #"when the data contains y-values in a column, use stat="identity"
  geom_errorbar(aes(ymax = mnormavg + sdnormavg, ymin=mnormavg - sdnormavg), width=0.1)+
  theme(panel.background = element_rect(fill="white"))+
  theme(axis.ticks = element_line(color="black", size = 0.2))+
  theme(axis.text.y = element_text(color="black", size = 10))+
  theme(axis.text.x = element_text(color="black", size = 14))+
  theme(axis.line = element_line(color="black", size = 0.2))+
  theme(axis.title.y = element_text(size = 16))+
  ylab("Normalized Luceriferase Output")+
  xlab(NULL)

#ordered factors in the order that I ultimately want them on the geom_point() plot below.
normdat.dt.noaug21may18[ ,gene_f:=factor(gene, levels=c("mak2", "nrc1", "ham5", "bem1", "pp1", "ada3", "adv1", 
                                                        "so", "ham4", "ham8", "ham9", "ham7", "ham6", "nox1", 
                                                        "ham11", "ham14", "mob3", "spr7", "ptp2", "vib1", 
                                                        "nik2", "cse1", "lao1", "arg15", "ham12", "WT"))]
#Discrete scatterplot of total normalized Luciferase data with aug21 and may18 experiments:
normdat.dt.m <-melt(normdat.dt.noaug21may18, 
                    id.vars=c("date", "strain", "gene", "isolate", "group", "group_f", "gene_f"), 
                    na.rm=TRUE)
ggplot(normdat.dt.m)+
  geom_point(aes(gene, value, color=group_f))+
  scale_x_discrete(limits=normdat.dt.m[order(gene_f), unique(gene)])+
  #scale_x_discrete(limits=normdat.dt.m[order(gene), unique(gene)])+ #alphabetical by gene
  theme(axis.text.x = element_text(angle=90, color="black"))+
  theme(legend.title=element_blank())+
  theme(panel.background = element_rect(fill="white"))+
  theme(panel.grid.minor = element_blank())+
  theme(panel.grid.major = element_line(color="grey90"))+  #grey100 = white, grey0 = black
  theme(axis.ticks = element_line(color="grey90"))+
  theme(axis.text.y = element_text(color="black"))+
  theme(legend.key = element_rect(fill = "white", color = "white"))+
  theme(legend.background = element_rect(colour = "grey70",size=0.25))+
  ylab("Normalized Luceriferase Output")+
  xlab(expression("Gene"))

#Discrete scatterplot of total normalized data, without any experiments removed.
normdat.dt[ ,group_f:=factor(group, levels=c('low', 'med', 'high'))]
normdat.dt[ ,gene_f:=factor(gene, levels=c("mak2", "nrc1", "ham5", "bem1", "pp1", "ada3", "adv1", 
                                           "so", "ham9", "ham8", "ham4", "ham7", "ham6", "nox1", 
                                           "ham11", "ham14", "mob3", "ras2", "pp2A", "mik1", 
                                           "amph1", "spr7", "ptp2", "vib1", "nik2", "cse1", 
                                           "lao1", "arg15", "doc12", "ham12", "WT"))]
normdat.dt.m <-melt(normdat.dt, 
                    id.vars=c("date", "strain", "gene", "isolate", "group", "group_f", "gene_f"), 
                    na.rm=TRUE)
ggplot(normdat.dt.m)+
  geom_point(aes(gene, value, color=group_f))+
  scale_x_discrete(limits=normdat.dt.m[order(gene_f), unique(gene)])+
  #scale_x_discrete(limits=normdat.dt.m[order(gene), unique(gene)])+ #alphabetical by gene
  theme(axis.text.x = element_text(angle=90, color="black"))+
  theme(legend.title=element_blank())+
  theme(panel.background = element_rect(fill="white"))+
  theme(panel.grid.minor = element_blank())+
  theme(panel.grid.major = element_line(color="grey90"))+  #grey100 = white, grey0 = black
  theme(axis.ticks = element_line(color="grey90"))+
  theme(axis.text.y = element_text(color="black"))+
  theme(legend.key = element_rect(fill = "white", color = "white"))+
  theme(legend.background = element_rect(colour = "grey70",size=0.25))+
  ylab("Normalized Luceriferase Output")+
  xlab(expression("Gene"))
  


#Sidenote: I can use the path instead of file.choose() to import data,

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
