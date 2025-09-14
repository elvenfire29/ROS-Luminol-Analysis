#Version 2025. September 04.
#ROS (Luminol) Assay Analysis:
#Averages per time point, Sums of Treatments for all Time Points, ANOVA comparing Sums


#PLEASE READ
#Need packages: ggplot2  dplyr  magrittr  stringr  filesstrings


#Pre-Input Instructions:
#1. Download the above packages into RStudio if you have not previously done so.
#2. Prep input by deleting all cells other than the fluorescence data with the Time row and accompanying column.
#3. Ensure there are no skipped or empty wells. Every well for the 96-well plate at every time point should have a value.
    #Please note that zero is a valid value and will not affect the calculations as it will be skipped when other data for the treatment is available.
#4. Save as CSV file.


########## User Input Below ############

#User Input I: Steps 1-3 below must be changed to tailor the analysis to a specific data set on a local computer.

#Step 1: Select your working directory (the folder with the input file)
setwd("C:/Users/pippi/Downloads")


#Step 2: Select your input CSV file
plate_no <- 'ROS_flg22.csv'


#Step 3: Name the treatments and/or lines used.
#Per specific experiment, please change sample names accordingly from A to H, OR from 1 to 12 as on the plate.
#The number and order of labels matters!! This is because the code counts the number (needs 8 or 12) and uses the order.
#Please do NOT use duplicate treatment names. Does not combine to a single group.
#Please avoid using a backslash in any treatment names!! Please also avoid "NA" as a label when using ANOVA!
treatment <- c(
  'Sample1',	'Sample2', 'Sample3',	'Sample4',	'Sample5',
  'Sample6', 'Sample7', 'Sample8'
)


#################
#User Input II: Options 1-4 can be changed per user's specific needs.

#Option 1: Run a One-Way ANOVA test with Tukey HSD post-hoc?
#Code will run the test if TRUE is selected. If FALSE is selected, it will not.
ANOVA <- TRUE

#################
#Option 2: Run a two-sided t-test?
#If yes, type TRUE; if no, type FALSE.
#If an error message says data is too similar, write FALSE and use ANOVA or do t-test by hand.
t_test <- TRUE

#################
#Option 3: Include graphs?
#If you want graphs of the fluorescence curves and a bar plot comparing the total sum of fluorescence, type TRUE. Otherwise, type FALSE.
GRAPHS <- TRUE

#Option 3 addition: Graph with SEM or SD?
#If you want to use SEM (standard error of mean) on the bar plots, type TRUE
#If you want to use SD (standard deviation) on the bar plots, type FALSE
SEM_CHOOSE <- TRUE

#################
#Option 4: How does your input file look?
#If your plate reader output lists wells by A1,A2,A3.... then type TRUE.
#If your reader output lists A1,B1,C1.... then type FALSE.
orientation <- FALSE

#################



############################################
#END OF USER INPUT
########################################################################################################
########################################################################################################
########################################################################################################
#############################################################################

########################
#Call libraries#
#######libraries#################
library(ggplot2)
library(dplyr)
library(stringr)
library(filesstrings)
library(magrittr)

#########################
# FUNCTION DEFINITIONS #
##########functions##############

GetTimepoints <- function(names) {
  #get actual timepoints
  timepoints <- c()

  #testing timepoints for minutes by splitting at h and looking at remaining lengths
  for (n in names) {
    header <- n[1]
    time <- 0

    #if it starts with an "X", remove it
    if (substr(header, 1, 1) == "X") {
      header <- substr(header, 2, nchar(header))
    }

    #if it end in " min", remove that whole chunk, and add the minutes back
    if (substr(header, nchar(header) - 2, nchar(header)) == "min") {
      #strip the "min" part
      header <- substr(header, 0, nchar(header) - 4)
      #split on the "h" if it exists
      if (substr(header,nchar(header)-1,nchar(header))=="h"){
      sepTime <- unlist(strsplit(header, ".h."))
      header <- (as.numeric(sepTime[1])*60)
      time <- time + as.numeric(sepTime[2])
      }
    }

    #if it ends in " h", remove that
    if (substr(header, nchar(header), nchar(header)) == "h") {
      header <- substr(header, 0, nchar(header) - 2)
      header <- (as.numeric(header[1])*60)
    }

    #we should be left with just the time in minutes, so let's add them
    time <- time + as.numeric(header)
    time <- na.omit(time)

    #append it to the timepoints array
    timepoints <- cbind(timepoints, time)
  }

  return(timepoints)
}


GetWellSeq <- function(i, treatment, orientation) {
  #change collection of samples based on the plate readout
  if (orientation == TRUE) {
    seq_12 <- seq(i, 84 + i, 12)
    seq_8 <- seq((i - 1) * 12 + 1, 12 * i, 1)
  }
  #change collection of samples based on the plate readout
  if (orientation == FALSE) {
    seq_12 <- seq(i, 88 + i, 8)
    seq_8 <- seq((i - 1) * 8 + 1, 8 * i, 1)
  }
  if (length(treatment) == 8) {
    return(seq_8)
  } else {
    return(seq_12)
  }
}

################################code starts/functions end#################################################

#get working directory
wd <- c(getwd())

#gets labels (time points) from infile
time_points = read.csv(as.character(plate_no))
time_points = time_points[, 2:ncol(time_points)]
names = colnames(time_points)
names = as.matrix(names)

#Get times in minutes
ZTtim<-GetTimepoints(names)

#not 8 or 12 treatments. doesn't know how to handle that
if ((length(treatment) != 8) & (length(treatment) != 12)) {
  stop("Please input labels for 8 or 12 treatments.")
}

#program decides if there are 12 or 8 samples per treatment
if (length(treatment) == 12) {
  print("There are 12 treatments.")
} else {
  print("There are 8 treatments.")
}

#repetitions per treatment changes based on number of treatments
rep <- if (length(treatment) == 8){12}else{8}

plate <- plate_no


#make treatment labels ok for file names by taking out certain symbols /  :  -  .  ~
#the backslash is never ok to use here

#Save the original treatment labels for use inside of files
treatment0<-treatment

#Remove the trouble causing symbols and replace them with a safe symbol
for (i in 1:length(treatment)){
  treatment[i]<-gsub('/', '_', treatment[i], ignore.case = FALSE, perl = FALSE, fixed = TRUE, useBytes = FALSE)
  treatment[i]<-gsub(':', '_', treatment[i], ignore.case = FALSE, perl = FALSE, fixed = TRUE, useBytes = FALSE)
  treatment[i]<-gsub('-', '_', treatment[i], ignore.case = FALSE, perl = FALSE, fixed = TRUE, useBytes = FALSE)
  treatment[i]<-gsub('.', '_', treatment[i], ignore.case = FALSE, perl = FALSE, fixed = TRUE, useBytes = FALSE)
  treatment[i]<-gsub('~', '_', treatment[i], ignore.case = FALSE, perl = FALSE, fixed = TRUE, useBytes = FALSE)
}

#remove repeated treatment names
treatment_dedup_2 <- treatment[!duplicated(treatment)]
treatment_dedup <- treatment0[!duplicated(treatment0)]

###write dividers for eventual outfiles
do <- vector()
for (i in 1:length(time_points)) {
  du <- c(".")
  do <- cbind(do, du)
}
Average <- as.vector(do)
SEM <- as.vector(do)

da <- vector()

for (i in 1:length(treatment0)){
  de <- c(".")
  da<-cbind(da,de)
}
IndividualWellSums <- as.vector(da)
TreatmentTotalSums <- as.vector(da)
SEMofTreatmentSums <- as.vector(da)
SEMofIndividualWells <- as.vector(da)

##############values############
#get unique treatment names
treat_unique<-unique(treatment0)

#dataset reorganizing so that treatments and time points are grouped for calculations
dat1 <- read.csv(plate)
xnum <- dat1[, 1]
#taking out "sample x" column
dat1 <- dat1[, -1]

#change into individual time points for each treatment (this is specific to treatment number)
breo <- c(seq(1:rep))
for (i in 1:length(treatment0)) {
  bery <- vector()
  for (x in 1:length(names)) {
    bver<-vector()
    well_seq <- GetWellSeq(i,treatment0,orientation)
    bver <- unlist(dat1[well_seq, x])
    if (x == 1) {
      bery <- c(bver)
    }
    if (x != 1) {
      bery <- cbind(bery, bver)
    }
  }
  breo <- cbind(breo, bery)
}
breo <- breo[, -1]

#average each column (data for each time point) or collect it for individual well time points
#set up vectors
averageTP <- vector()
combined_sd <- vector()
indiv <- vector()
combinedTreatmentTP <- vector()

for (x in 1:length(colnames(breo))) {

  #get rid of NA values
  temp3 <- breo[, x]
  temp4 <- na.omit(temp3)

  temp<-vector()
  temp<-c(temp4)

  #get rid of 0 values unless 0 is all there is
  did<-cbind.data.frame(temp,temp)
  if (sum(did[])!=0){
  ans = did[rowSums(did[])>0,]
  temp<-ans[,1]
  }

  #collect individual well values
  indiv <- cbind(indiv,temp4)

  #do calculations
  #find average (mean)
  av_time <- mean(temp)
  averageTP <- cbind(averageTP, av_time)

  #do calculation for standard error of mean (SEM)
  sam_treatment <- (sd(temp) / sqrt(length(temp)))
  combined_sd <- cbind(combined_sd, sam_treatment)
}

#find the total sum and SEM across time for each replicate
#Set vectors
temp11<-vector()
wellSums<-vector()
SumSEMI <-vector()
Individuals<-vector()

for (i in 1:length(treatment0)){
  #get rid of any NA values
  temp11 <- na.omit(breo)

  temp12<-vector()
  temp14<-vector()
  tempI<-vector()

  for (x in 1:rep){
    if (i==1){
      #calculations
      #sum
      sums <- breo[x,seq((i-1)*length(names)+1,(i-1)*length(names)+length(names),1)]
      #get rid of 0 values unless 0 is all there is
      did<-cbind.data.frame(sums,sums)
      if (sum(did[])!=0){
        ans = did[rowSums(did[])>0,]
        sums<-ans[,1]
      }
      sums<-sum(sums)
      temp12 <- rbind(temp12,sums)

      #sem
      semsI <-breo[x,seq((i-1)*length(names)+1,(i-1)*length(names)+length(names),1)]
      #get rid of 0 values unless 0 is all there is
      did<-cbind.data.frame(semsI,semsI)
      if (sum(did[])!=0){
        ans = did[rowSums(did[])>0,]
        semsI<-ans[,1]
      }
      semsI <- sd(semsI)/sqrt(length(semsI))
      temp14<-rbind(temp14,semsI)

      #indiv
      tempV<-breo[x,seq((i-1)*length(names)+1,(i-1)*length(names)+length(names),1)]
      #tempV<-c(indiv[x,seq(i,ncol(indiv)-((1-(i/length(treatment)))*length(treatment)),length(treatment))])
      tempI<-rbind(tempI,tempV)
    }
    else{
      #calculations
      #sum
      sums <- breo[x,seq((i-1)*length(names)+1,(i-1)*length(names)+length(names),1)]
      #get rid of 0 values unless 0 is all there is
      did<-cbind.data.frame(sums, sums)
      if (sum(did[])!=0){
        ans = did[rowSums(did[])>0,]
        sums<-ans[,1]
      }
      sums <-sum(sums)
      temp12 <- rbind(temp12,sums)

      #sem
      semsI <-breo[x,seq((i-1)*length(names)+1,(i-1)*length(names)+length(names),1)]
               #get rid of 0 values unless 0 is all there is
               did<-cbind.data.frame(semsI, semsI)
               if (sum(did[])!=0){
                 ans = did[rowSums(did[])>0,]
                 semsI<-ans[,1]
               }
      semsI<- sd(semsI)/sqrt(length(semsI))
      temp14<-rbind(temp14,semsI)

      #indiv
      tempV<-breo[x,seq((i-1)*length(names)+1,(i-1)*length(names)+length(names),1)]
      #tempV<-c(indiv[x,seq(i,ncol(indiv)-((1-(i/length(treatment)))*length(treatment)),length(treatment))])
      tempI<-rbind(tempI,tempV)
    }
  }

  if(i==1){
    wellSums <- c(temp12)
    SumSEMI <- c(temp14)
    #individual wells
    Individuals<-rbind(Individuals,tempI)
  }
  else{
    wellSums <- cbind(wellSums,temp12)
    SumSEMI <- cbind(SumSEMI,temp14)
    Individuals<-rbind(Individuals,tempI)
  }
}

#label indiv well tps
colnames(Individuals)<-ZTtim
#get labels for replicates
Treatments <-vector()
for (i in treatment0){
  for (x in 1:rep){
    Treatments<-rbind(Treatments,i)
  }
}
Individuals <- cbind(treatment0,Individuals)
rownames(Individuals)<-NULL

#move labels if same-name not consecutive
reorgNames<-vector()
if (length(treat_unique)!=length(treatment0)){
for (i in 1:length(treatment0)){
  for (n in 1:length(Treatments)){
    #grab all treatments with same name
    if(str_detect(Treatments[n], treatment0[i])){
      reorgNames<-rbind(reorgNames,Individuals[n,])
    }
  }
}
  Individuals<-c(reorgNames)
  rownames(Individuals)<-NULL
}

#reorganize averages and SEM again so that the values for each treatment at each time point go top to bottom
graphData <- vector()
SEMdata <- vector()

for (i in 1:length(treatment0)) {
  reorgan <-
    unlist(averageTP[seq((i - 1) * length(names) + 1,
                         (i - 1) * length(names) + length(names),
                         1)])
  reorganSEM <-
    unlist(combined_sd[seq((i - 1) * length(names) + 1,
                           (i - 1) * length(names) + length(names),
                           1)])

  graphData <- rbind(graphData, reorgan)
  SEMdata <- rbind(SEMdata, reorganSEM)
}

#make duplicates so we don't alter the originals that may be used for graphing
averages <- graphData
SEMvalues <- SEMdata

#Find the total sum and SEM across time for each treatment
#uses the averages of each well to find SEM and final sum
#set up vectors and remove NA values
temp16 <- na.omit(averages)
TreatmentSums <-vector()
SumSEM <- vector()
SumSD <- vector()

for (i in 1:length(treatment)){
  #sums
  Totalsums<-sum(temp16[i,])
  TreatmentSums <- cbind(TreatmentSums,Totalsums)

  #sem
  sems <-(sd(temp16[i,]) / sqrt(length(temp16[i,])))
  SumSEM <- cbind(SumSEM,sems)
  SumSD <- cbind(SumSD,sd(temp16[i,]))
}

#label rows
row.names(averages) <- treatment0
row.names(SEMvalues) <- treatment0
row.names(wellSums) <- seq(1,rep,1)
row.names(SumSEMI) <- seq(1,rep,1)

#label columns
colnames(averages) <- ZTtim
colnames(SEMvalues) <- ZTtim
colnames(wellSums) <- treatment0
colnames(TreatmentSums) <- treatment0
colnames(SumSEMI) <- treatment0

#############most csv outfiles and outfile folder#############
#create new folder for output
plate_num <- strsplit(plate_no, split = ".csv")
folder <- paste(plate_num[1], 'output')
if (!dir.exists(folder)) {
  dir.create(folder)
}

#set location of new sub-folders
wd <- getwd()
newloco <- paste(wd, folder, sep = "/")

#create sub folders for output
folder2 <- paste(plate_num[1], 't_test')
folder3 <- paste(plate_num[1], 'ANOVA')
folder4 <- paste(plate_num[1], 'PRISM and Graphing Data')

t_testFolder <- paste(newloco,folder2,sep='/')
ANOVAFolder <- paste(newloco,folder3,sep='/')
PrismFolder <- paste(newloco,folder4,sep='/')

#make ANOVA directory
if(ANOVA){
  if (!dir.exists(ANOVAFolder)) {
    dir.create(ANOVAFolder)
  }
}

if(t_test){
if (!dir.exists(t_testFolder)) {
  dir.create(t_testFolder)
}
}

if (!dir.exists(PrismFolder)) {
  dir.create(PrismFolder)
}

################making some out-files##############################

#add total ROS counts per seedling and per genotype or treatment
finalCount <- rbind(TreatmentTotalSums, TreatmentSums, SEMofTreatmentSums, SumSEM, IndividualWellSums, wellSums,SEMofIndividualWells,SumSEMI)
otherout <- paste(plate_num[1], 'Total ROS Counts.csv')
write.csv(x = finalCount, otherout)
#move files to new folder
file.move(otherout, newloco, overwrite = TRUE)

#wells averaged stats outfile
  final <- rbind(Average, averages, SEM, SEMvalues)
  #set up name of outfile
  out_file <- paste(plate_num[1], 'Averages Over Time.csv')
  #write outfile
  write.csv(x = final, out_file)
  file.move(out_file, PrismFolder, overwrite = TRUE)

  #individual wells outfile
  out<-paste(plate_num[1],"Individuals Over Time.csv")
  #write outfile
   write.csv(x = Individuals, out)
   file.move(out, PrismFolder, overwrite = TRUE)

  #transpose averages
  PRISM_DAT2<-t(averages)
  file2<-paste(plate_num[1],"PRISM Ready Averages.csv")
  write.csv(x = PRISM_DAT2, file2)
  file.move(file2, PrismFolder, overwrite = TRUE)

  #PRISM ready data for sums
  file3<-paste(plate_num[1],"PRISM Ready Sums.csv")
  write.csv(x = wellSums, file3)
  file.move(file3, PrismFolder, overwrite = TRUE)

  #################################
  #make sure treatment names don't duplicate in data used for individual graphs
  treatment_dedup <- treatment0[!duplicated(treatment0)]
  #labels for data
  rownames(graphData) <- treatment_dedup
  colnames(graphData) <- ZTtim


  ###################GRAPHS##############################
#outputs graphs as PDF
 if (GRAPHS){
  #make file name
  graph <- paste(plate_num[1], 'Graphs.pdf')

  #########################################start pdf#####
  #Open a pdf file
  pdf(graph, width = 9, height = 35)

  #2 columns and TBD rows of graphs
  par(
    mfrow = c(ceiling(length(row.names(graphData)) + 2), 2),
    las = 2,
    #set outer (outside graph block) margins: below, left, top, right
    oma = c(1, 1, 1, 1),
    #set inner (inter-graph) margins: below, left, top, right
    mar=c(9,11,2,2)
  )

  #define mround
  mround <- function(number, multiple) multiple * ceiling(number/multiple)

  #################################Graph of Sums###############################
   #get upper y limit for axis
    maxP <- max(TreatmentSums)
    lim <- c(0)
    if (maxP <= 500) {
      lim <- mround(maxP, 100)
    }
    if (maxP > 500) {
      lim <- mround(maxP, 500)
    }
    if (maxP >= 1000) {
      lim <- mround(maxP, 1000)
    }
    if (maxP >= 5000) {
      lim <- mround(maxP, 3000)
    }
    if (maxP >= 10000) {
      lim <- mround(maxP, 7000)
    }
    if (maxP >= 15000) {
      lim <- mround(maxP, 10000)
    }
    if (maxP >= 20000) {
      lim <- mround(maxP, 10000)
    }
    if (maxP >= 35000) {
      lim <- mround(maxP, 15000)
    }
    if (maxP >= 50000) {
      lim <- mround(maxP, 20000)
    }
    if (maxP >= 75000) {
      lim <- mround(maxP, 25000)
    }
    if (maxP >= 100000) {
      lim <- mround(maxP, 30000)
    }
    if (maxP >= 125000) {
      lim <- mround(maxP, 35000)
    }
    if (maxP >= 150000) {
      lim <- mround(maxP, 50000)
    }
    if (maxP >= 200000) {
      lim <- mround(maxP, 55000)
    }
    if (maxP >= 250000) {
      lim <- mround(maxP, 80000)
    }
    if (maxP >= 300000) {
      lim <- mround(maxP, 100000)
    }
    if (maxP >= 350000) {
      lim <- mround(maxP, 150000)
    }
    if (maxP >= 400000) {
      lim <- mround(maxP, 200000)
    }
    if (maxP >= 450000) {
      lim <- mround(maxP, 300000)
    }
    if (maxP >= 500000) {
      lim <- mround(maxP, 350000)
    }
    if (maxP >= 550000) {
      lim <- mround(maxP, 400000)
    }

    #get title for graph
    title2 <- c("ROS Sums")
    par(mgp=c(5,1,0))

    #try to avoid scientific notation on the y axis
    options(scipen=10)

    #graph
    barplot(
      height = TreatmentSums,
      names = colnames(TreatmentSums),
      xlab = NULL,
      ylab = "Fluorescence Units",
      ylim=c(0,lim),
      main = title2,
      horiz=FALSE,
      axes=TRUE,
      mar=c(10,10,2,2),
      las=2,
      col = "purple"
    )

    #add a line on the top of each bar for easier comparison
    #abline(h=TreatmentSums,col="plum")

    #Graph it again, but now with error bars
    ROS2<-barplot(
      height = TreatmentSums,
      names = colnames(TreatmentSums),
      xlab = NULL,
      ylab = "Fluorescence Units",
      ylim=c(0,lim),
      main = paste(title2,"with Error Bars"),
      horiz=FALSE,
      axes=TRUE,
      las=2,
      col = "purple"
    )
    #get the error bars
    if(SEM_CHOOSE==TRUE){
    #use SEM
    #error bars hack: we draw arrows but with very special "arrowheads"
    arrows(
      ROS2,
      as.numeric(TreatmentSums) - as.numeric(SumSEM),
      ROS2,
      as.numeric(TreatmentSums) + as.numeric(SumSEM),
      length = 0.05,
      angle = 90,
      code = 3,
      col = 'black'
    )
    }else{
      #use SD instead of SEM
      #error bars hack: we draw arrows but with very special "arrowheads"
      arrows(
        ROS2,
        as.numeric(TreatmentSums) - as.numeric(SumSD),
        ROS2,
        as.numeric(TreatmentSums) + as.numeric(SumSD),
        length = 0.05,
        angle = 90,
        code = 3,
        col = 'black'
      )
    }

   #############################Graphs of each Treatment################
  for (i in 1:length(row.names(graphData)))
  {
    #get upper y limit for axis
    maxP <- max(graphData[i,])
    lim <- c(0)
    if (maxP <= 500) {
      lim <- mround(maxP, 100)
    }
    if (maxP > 500) {
      lim <- mround(maxP, 500)
    }
    if (maxP >= 1000) {
      lim <- mround(maxP, 1000)
    }
    if (maxP >= 5000) {
      lim <- mround(maxP, 3000)
    }
    if (maxP >= 10000) {
      lim <- mround(maxP, 7000)
    }
    if (maxP >= 15000) {
      lim <- mround(maxP, 10000)
    }
    if (maxP >= 20000) {
      lim <- mround(maxP, 15000)
    }
    if (maxP >= 35000) {
      lim <- mround(maxP, 15000)
    }
    if (maxP >= 50000) {
      lim <- mround(maxP, 20000)
    }
    if (maxP >= 75000) {
      lim <- mround(maxP, 25000)
    }
    if (maxP >= 100000) {
      lim <- mround(maxP, 30000)
    }
    if (maxP >= 125000) {
      lim <- mround(maxP, 35000)
    }
    if (maxP >= 150000) {
      lim <- mround(maxP, 40000)
    }
    if (maxP >= 200000) {
      lim <- mround(maxP, 50000)
    }
    if (maxP >= 250000) {
      lim <- mround(maxP, 80000)
    }
    if (maxP >= 300000) {
      lim <- mround(maxP, 100000)
    }
    if (maxP >= 350000) {
      lim <- mround(maxP, 150000)
    }
    if (maxP >= 400000) {
      lim <- mround(maxP, 200000)
    }
    if (maxP >= 450000) {
      lim <- mround(maxP, 300000)
    }
    if (maxP >= 500000) {
      lim <- mround(maxP, 350000)
    }
    if (maxP >= 550000) {
      lim <- mround(maxP, 400000)
    }

    #get title for graph
    title <- paste(treatment_dedup[i], "ROS")

    #graph
    plot(
      x = ZTtim,
      y = graphData[i,],
      ylim = c(0, lim),
      xlim = c(0, length(names) + 5),
      xlab = "Time (minutes)",
      ylab = "Fluorescence Units",
      main = title,
      type = 'p',
      col = "blue",
      xaxt = 'n'
    )
    lines(ZTtim, graphData[i,], type = 'o', col = 'blue')
    #x axis
    if (length(names) < 100) {
      axis(1, c(seq(12, length(names) + 5, 12)))
    }
    if (length(names) >= 100 &
        length(names) < 200) {
      axis(1, c(seq(12, length(names) + 5, 24)))
    }
    if (length(names) >= 200) {
      axis(1, c(seq(24, length(names) + 5, 48)))
    }

    #y axis
if (lim !=0){
  if (lim>=200000){
    axis(2, c(seq(70000, mround(lim, 70000), 70000)))
  }
  else if (lim>=100000){
    axis(2, c(seq(50000, mround(lim, 50000), 50000)))
  }
      else if (lim>=60000){
        axis(2, c(seq(20000, mround(lim, 20000), 20000)))
      }
      else if (lim >= 45000) {
        axis(2, c(seq(10000, mround(lim, 10000), 15000)))
      }
      else if (lim >= 20000) {
        axis(2, c(seq(5000, mround(lim, 5000), 10000)))
      }
      else if (lim >= 10000) {
        axis(2, c(seq(2000, mround(lim, 2000), 5000)))
      }
      else if (lim >= 5000) {
        axis(2, c(seq(2000, mround(lim, 1000), 2000)))
      }
      else if (lim >= 2500) {
        axis(2, c(seq(500, mround(lim, 500), 1000)))
      }
      if (2500 > lim) {
        axis(2, c(seq(250, mround(lim, 250), 250)))
      }
      if (2500 >= lim) {
        axis(2, c(seq(250, mround(lim, 250), 250)))
      }
      if (250 >= lim) {
        axis(2, c(seq(50, mround(lim, 50), 50)))
      }
    if (100 >= lim) {
      axis(2, c(seq(20, mround(lim, 10), 20)))
    }
    if (50 >= lim) {
      axis(2, c(seq(10, mround(lim, 5), 10)))
    }
    if (20 >= lim) {
      axis(2, c(seq(5, mround(lim, 5), 5)))
    }
}
    if (lim == 0){
      axis(2,c(seq(0, 10, 1)))
    }

  #Do it again, but now with error bars
    #get the error bars
    well_seq <- GetWellSeq(i, treatment, orientation)
    error_bar_points <- dat1[well_seq,]
    error_bars <- c()
    for (j in 1:ncol(error_bar_points)) {
      error_bar_cur_column <- error_bar_points[[i]]
      error_bar_value <-
        sd(error_bar_cur_column) / sqrt(length(error_bar_cur_column))
      error_bars <- append(error_bars, error_bar_value)
    }

    #get upper y limit for axis
    maxP <- max(graphData[i,])
    lim <- c(0)
    if (maxP <= 500) {
      lim <- mround(maxP, 100)
    }
    if (maxP > 500) {
      lim <- mround(maxP, 500)
    }
    if (maxP >= 1000) {
      lim <- mround(maxP, 1000)
    }
    if (maxP >= 10000) {
      lim <- mround(maxP, 5000)
    }
    if (maxP >= 20000) {
      lim <- mround(maxP, 10000)
    }

    #get title for graph
    title <- paste(treatment_dedup[i], "ROS")

    #graph
    plot(
      x = ZTtim,
      y = graphData[i,],
      xlim = c(0, length(names) + 5),
      ylim = c(0, lim),
      xlab = "Time (minutes)",
      ylab = "Fluorescence Units",
      main = title,
      type = 'p',
      col = "blue",
      xaxt = 'n'
    )

    #add arrow bars [sic]
    arrows(
      ZTtim,
      graphData[i,] - error_bars,
      ZTtim,
      graphData[i,] + error_bars,
      length = 0.05,
      angle = 90,
      code = 3
    )

    lines(ZTtim, graphData[i,], type = 'o', col = 'blue')
    #x axis
    if (length(names) < 100) {
      axis(1, c(seq(12, length(names) + 5, 12)))
    }
    if (length(names) >= 100 &
        length(names) < 200) {
      axis(1, c(seq(12, length(names) + 5, 24)))
    }
    if (length(names) >= 200) {
      axis(1, c(seq(24, length(names) + 5, 48)))
    }
    #y axis
    if (lim !=0){
      if (lim>=200000){
        axis(2, c(seq(70000, mround(lim, 70000), 70000)))
      }
      else if (lim>=100000){
        axis(2, c(seq(50000, mround(lim, 50000), 50000)))
      }
      else if (lim>=60000){
        axis(2, c(seq(20000, mround(lim, 20000), 20000)))
      }
      else if (lim >= 45000) {
        axis(2, c(seq(10000, mround(lim, 10000), 15000)))
      }
      else if (lim >= 20000) {
        axis(2, c(seq(5000, mround(lim, 5000), 10000)))
      }
      else if (lim >= 10000) {
        axis(2, c(seq(2000, mround(lim, 2000), 5000)))
      }
      else if (lim >= 5000) {
        axis(2, c(seq(2000, mround(lim, 1000), 2000)))
      }
      else if (lim >= 2500) {
        axis(2, c(seq(500, mround(lim, 500), 1000)))
      }
      if (2500 > lim) {
        axis(2, c(seq(250, mround(lim, 250), 250)))
      }
      if (2500 >= lim) {
        axis(2, c(seq(250, mround(lim, 250), 250)))
      }
      if (250 >= lim) {
        axis(2, c(seq(50, mround(lim, 50), 50)))
      }
      if (100 >= lim) {
        axis(2, c(seq(20, mround(lim, 10), 20)))
      }
      if (50 >= lim) {
        axis(2, c(seq(10, mround(lim, 5), 10)))
      }
      if (20 >= lim) {
        axis(2, c(seq(5, mround(lim, 5), 5)))
      }
    }
    if (lim == 0){
      axis(2,c(seq(0, 10, 1)))
    }
  }

  ###now make a graph with all the lines
  maxP <- max(graphData)
  lim <- c(0)
  if (maxP <= 500) {
    lim <- mround(maxP, 100)
  }
  if (maxP > 500) {
    lim <- mround(maxP, 500)
  }
  if (maxP >= 1000) {
    lim <- mround(maxP, 1000)
  }
  if (maxP >= 10000) {
    lim <- mround(maxP, 5000)
  }
  if (maxP >= 20000) {
    lim <- mround(maxP, 10000)
  }

  title <- paste("ROS in", plate_num[1])
  plot(
    x = ZTtim,
    y = graphData[1,],
    ylim = c(0, lim),
    xlim = c(0, length(names) + 5),
    xlab = "Time (minutes)",
    ylab = "Fluorescence Units",
    main = title,
    type = 'l',
    col = "blue",
    xaxt = 'n'
  )
  lines(ZTtim, graphData[1,], type = 'l', col = 'blue')

  #set axis scale
  #x axis
  if (length(names) < 100) {
    axis(1, c(seq(12, length(names) + 5, 12)))
  }
  if (length(names) >= 100 &
      length(names) < 200) {
    axis(1, c(seq(12, length(names) + 5, 24)))
  }
  if (length(names) >= 200) {
    axis(1, c(seq(24, length(names) + 5, 48)))
  }
  #y axis
  if (lim !=0){
    if (lim>=200000){
      axis(2, c(seq(70000, mround(lim, 70000), 70000)))
    }
    else if (lim>=100000){
      axis(2, c(seq(50000, mround(lim, 50000), 50000)))
    }
    else if (lim>=60000){
      axis(2, c(seq(20000, mround(lim, 20000), 20000)))
    }
    else if (lim >= 45000) {
      axis(2, c(seq(10000, mround(lim, 10000), 15000)))
    }
    else if (lim >= 20000) {
      axis(2, c(seq(5000, mround(lim, 5000), 10000)))
    }
    else if (lim >= 10000) {
      axis(2, c(seq(2000, mround(lim, 2000), 5000)))
    }
    else if (lim >= 5000) {
      axis(2, c(seq(2000, mround(lim, 1000), 2000)))
    }
    else if (lim >= 2500) {
      axis(2, c(seq(500, mround(lim, 500), 1000)))
    }
    if (2500 > lim) {
      axis(2, c(seq(250, mround(lim, 250), 250)))
    }
    if (2500 >= lim) {
      axis(2, c(seq(250, mround(lim, 250), 250)))
    }
    if (250 >= lim) {
      axis(2, c(seq(50, mround(lim, 50), 50)))
    }
    if (100 >= lim) {
      axis(2, c(seq(20, mround(lim, 10), 20)))
    }
    if (50 >= lim) {
      axis(2, c(seq(10, mround(lim, 5), 10)))
    }
    if (20 >= lim) {
      axis(2, c(seq(5, mround(lim, 5), 5)))
    }
  }
  if (lim == 0){
    axis(2,c(seq(0, 10, 1)))
  }


  rain <-
    c(
      'black',
      'cyan',
      'green2',
      'dark green',
      'red',
      'orange',
      'darkorchid',
      'yellow1',
      'mediumpurple',
      'red4',
      'thistle',
      'sienna3',
      'turquoise4',
      'gray33',
      'slateblue3',
      'darkgoldenrod',
      'deeppink'
    )
  colors <- c('blue')
  for (i in 2:length(row.names(graphData)))
  {
    #get upper y limit for axis
    maxP <- max(graphData[i,])
    lim <- c(0)
    if (maxP <= 500) {
      lim <- mround(maxP, 100)
    }
    if (maxP > 500) {
      lim <- mround(maxP, 500)
    }
    if (maxP >= 1000) {
      lim <- mround(maxP, 1000)
    }
    if (maxP >= 10000) {
      lim <- mround(maxP, 5000)
    }
    if (maxP >= 20000) {
      lim <- mround(maxP, 10000)
    }

    #graph
    lines(
      x = ZTtim,
      y = graphData[i,],
      ylim = c(0, lim),
      xlim = c(0, length(names) + 5),
      xlab = "Time (minutes)",
      ylab = "Fluorescence Units",
      main = title,
      type = 'l'
    )
    lines(ZTtim, graphData[i,], type = 'l', col = rain[i])
    color <- c(rain[i])
    colors <- cbind(colors, color)
  }

  #make legend to side
  plot(
    0,
    0,
    type = "n",
    bty = "n",
    xaxt = "n",
    yaxt = "n",
    xlab = "",
    ylab = ""
  )
  legend(
    'topleft',
    legend = c(treatment_dedup),
    col = c(colors),
    pch = 10,
    bty = 'n',
    ncol = 2
  )


  # Close the pdf file
  dev.off()

  ###move graph file into new folder###
  file.move(graph, newloco, overwrite = TRUE)
}


  ##########################t-tests##############################
  #does t test control vs treatment and extracts p-value
  if(t_test){
    #loop for all controls on the plate (run each treatment as a control)
    #collect and test the ROS sums of each treatment
    sver<-t(wellSums)

      for (z in 1:length(treatment)) {
        control20<-sver[z,]

        scollect <- vector()
        scollect2 <- vector()

        #compare control to all other treatments
         for (i in 1:length(treatment)) {
          temp20<-sver[i,]

    #unequal variance, alt=means are not the same (two sided)
    #Welch two sample t-test
    #p-value adj method: bonferroni
    p_valueS <- t.test(control20, temp20, var.equal = F, p.adjust.method="bonferroni",alternative = c("two.sided"))$p.value

    #collect all p-values for this control
    scollect <- rbind(scollect,p_valueS)


    #for pairwise, every treatment must be same-length
      data6<-c()

      #get data into the same space and differentiate it
    if(length(control20)==length(temp2)){
      data6<-rbind(control20,temp20)
      data6<-t(data6)
      colnames(data6)<-c("Sample1","Sample2")

      #pairwise comparisons using t-tests with paired groups (no pooled SD)
      #p-value adj method: bonferroni
      p_valueS2 <- pairwise.t.test(data6, c("Sample1","Sample2"), p.adjust.method="bonferroni", var.equal = F, paired=TRUE, alternative=c("two.sided"))$p.value
      #collect all p-values for this control
      scollect2 <- rbind(scollect2,p_valueS2)
      #Set flag to false
      Sorry <- FALSE
      }else{
          Sorry <- TRUE
    }
         }

        #set labels
        combined_p_valuesS <- "p-value"
        combined_p_valuesS2 <- "p-value"

        #add rowname labels to p-values
        scollect <- cbind(treatment0,scollect)

    #welch
        #label
        scollect <- cbind(combined_p_valuesS,scollect)

        #create outfiles
        s<-treatment
        number<-paste("#",z,"__",sep="")
        out_fileS <- paste(number,s[z],'Welch p_value.csv')

        #write outfile
        write.csv(x=scollect,out_fileS)
        #move to a subfolder for t-test results
        file.move(file=out_fileS, t_testFolder, overwrite=TRUE)


    #if we were able to compute the pairwise t-test, include them. Otherwise, don't.
       if(length(scollect2)>1){

         #label
      scollect2 <- cbind(treatment0,scollect2)
      scollect2 <- cbind(combined_p_valuesS2,scollect2)

      #create outfiles
      s2<-treatment
      number<-paste("#",z,"__",sep="")
      out_fileS2 <- paste(number,s2[z],'Pairwise p_value.csv')

      #write outfile
      write.csv(x=scollect2,out_fileS2)

      #move to a subfolder for t-test results
      file.move(file=out_fileS2, t_testFolder, overwrite=TRUE)
      }
      }

    #use flag
    if (Sorry){
      print("Sorry! Treatments are not of the same length and thus cannot be compared pairwise.")
    }
  }


  ##########################ANOVA statistical tests##############################
  #does ANOVA and extracts p-values
  if(ANOVA){

        #reorganize for ANOVA
        DataForANOVA<-vector()
        almostANOVA<-vector()
        treatmentX<-colnames(wellSums)

        hi<-vector()
        for (x in 1:length(treatment0)){
          wholeCol <- vector()
          hey<-vector()
          for (i in 1:rep){
            #check if treatment is only 0 values
            wholeCol <- sum(wellSums[,x])
            #get rid of 0 values if not
            if (wholeCol != 0){
            didi<-rbind(wellSums[i,x],wellSums[i,x])
              if (sum(didi)==0){
                hi<-hi+1
                hey<-hey+1
              }else{
                DataForANOVA<-rbind(DataForANOVA,wellSums[i,x])
                almostANOVA <- rbind(almostANOVA,treatmentX[x])
               }
            }else{
              DataForANOVA<-rbind(DataForANOVA,wellSums[i,x])
              almostANOVA <- rbind(almostANOVA,treatmentX[x])
            }
          }
        }

        #Label
        DataForANOVA<-cbind(almostANOVA,DataForANOVA)

        #reorganize the file
        #find number of unique treatments
        Treatments0<-unique(DataForANOVA[,1])
        Treatments<-DataForANOVA[,1]

        #move same treatment names so that they are near one another
        if(length(Treatments0) != length(treatment)){
          potkettle<-c()
          blacklist<-c()
          for(i in 1:length(Treatments)){
            #find number of replicates per treatment
            Repeats<-which(Treatments==Treatments[i])
            #get all reps for a single treatment
            if(i==1){
              count=0
              for(k in Repeats){
                if(count==0){
                  potkettle<-c(DataForANOVA[k,])
                  count=count+1
                }
                else{
                  potkettle<-rbind(potkettle,DataForANOVA[k,])
                }
              }
            }
            else{
              #make sure we aren't repeating treatment replicates
              if((Treatments[i] %in% blacklist)==FALSE){
                for(p in Repeats){
                  potkettle<-rbind(potkettle,DataForANOVA[p,])
                }
              }
            }
            blacklist<-c(blacklist,Treatments[i])
          }
          DataForANOVA<-potkettle
        }

        #Label the column names
        colnames(DataForANOVA)<-c("Plant_Line","ROS_Sum")
        row.names(DataForANOVA)<-NULL

        #alphabetize data by plant line
        DataForANOVA<-DataForANOVA[order(DataForANOVA[,'Plant_Line']),]

        #write file for running ANOVA code for easier re-running in case of changes or well deletions later
        plate_num <- strsplit(plate_no, split = ".csv")
        titleT<-paste("ANOVA_Data",".csv",sep="")
        title2<-paste(plate_num,titleT)
        write.csv(x=DataForANOVA,title2,quote=FALSE,row.names = FALSE)
        if(ANOVA==FALSE){
          file.move(title2, newloco, overwrite = TRUE)
        }

        ###Start ANOVA statistics
        #read in data
        my_data<-read.csv(title2,header=TRUE,colClasses=c("factor","numeric"))

        #check data is read in properly
        print(summary(my_data))

        #run the anova
        results <- aov(ROS_Sum~Plant_Line,data = my_data)

        #print out the summary of values from the anova
        print(summary(results))

        #ad-hoc test
        post_results<-TukeyHSD(results)

        #put results into outfile to read-in below
        otheroutA <- c('ANOVA_Tukey All Results.csv')
        stats<-c(post_results,"p-values of Sums")
        write.csv(x = stats, otheroutA)

        #read in the all results ANOVA file
        ANOVA_PVALUES<-read.csv(otheroutA)

        #edit so that it is easier to read (i.e. take out unnecessary columns and leave p-value)
        for (i in 1:length(ANOVA_PVALUES)){
          if(i==1){
            ANOVA_PVALUES<-ANOVA_PVALUES[,-seq(2,4,1)]
          }
          if(i!=1){
            ANOVA_PVALUES<-ANOVA_PVALUES[,-seq(i+1,i+4,1)]
          }
        }
        #Change column names to be more descriptive
        ANOVA_COL<-c('Treatments Compared','p-value of Sums')
        colnames(ANOVA_PVALUES)<-ANOVA_COL

        #Look for ANOVAs with the control treatment
        comparisons<-ANOVA_PVALUES[,1]
        #make sure that if we can tell wt from wt with chemical, etc
        #add spaces to ends of each label to attempt to do so
        for (n in 1:length(comparisons)){
          comparisons[n] <-paste(".",comparisons[n],".",sep="")
        }

        #Define unique treatment names
        OGTreat_unique<-unique(my_data[,1])

        #There are all the controls
        for(k in 1:length(OGTreat_unique)){
          cup<-c()
          #take out control and separating - so that only the treatment tested against control remains
          test1<-paste("-",OGTreat_unique[k],".",sep="")
          test2<-paste(".",OGTreat_unique[k],"-",sep="")

          #search for each treatment separately
          for(n in 1:length(comparisons)){
            #Some treatments may include the name of other treatments, but with added stuff
            chonks1 <-comparisons[n]
            temp<-c()

            #if TRUE, then add the row to the collection cup for the file
            if(str_detect(chonks1, test1)==TRUE){
              #make sure ID is 100% correct
              ID<-sub(test1,"",chonks1)
              #take out the "." that is left in ID at beginning or end
              #count occurrences of "." so that we don't ruin annotation of treatment
              COUNT<-sum(charToRaw(ID) == charToRaw('.'))
              if(COUNT==1){
                #remove "."
                tempID<-str_remove(ID,"[.]")
              }
              #remove only the first or last "."
              if(COUNT>1){
                #find if "." is first or last
                if(unlist(gregexpr('.', ID))[1]==1){
                  #get rid of first "." only
                  tempID<-sub('.', '', ID)
                }
                last<-tail(unlist(gregexpr('.', ID)), n=1)
                temptemptemp<-unlist(ID)
                if(nchar(temptemptemp)==last){
                  #get rid of last "." only by knocking out the last character, which we assume is "."
                  tempID<-str_sub(ID,1,-2)
                }
              }
              #compare to other names to see if same
              for (u in 1:length(OGTreat_unique)){
                if(tempID==OGTreat_unique[u]){
                  if(tempID!=OGTreat_unique[k]){
                    temp<-ANOVA_PVALUES[n,]
                    if(length(cup)==0){
                      cup<-temp
                    }
                    else{
                      cup<-rbind(cup,temp)
                    }
                  }
                }
              }
            }
            #if the other way is TRUE, add it in!
            if(str_detect(chonks1, test2)==TRUE){
              #make sure ID is 100% correct
              ID<-sub(test2,"",chonks1)
              #take out the "." that is left in ID at beginning or end
              #count occurrences of "." so that we don't ruin annotation of treatment
              COUNT<-sum(charToRaw(ID) == charToRaw('.'))
              if(COUNT==1){
                #remove "."
                tempID<-str_remove(ID,"[.]")
              }
              #remove only the first or last "."
              if(COUNT>1){
                #find if "." is first or last
                if(unlist(gregexpr('.', ID))[1]==1){
                  #get rid of first "." only
                  tempID<-sub('.', '', ID)
                }
                last<-tail(unlist(gregexpr('.', ID)), n=1)
                temptemptemp<-unlist(ID)
                if(nchar(temptemptemp)==last){
                  #get rid of last "." only by knocking out the last character, which we assume is "."
                  tempID<-str_sub(ID,1,-2)
                }
              }
              #compare to other names to see if same
              for (u in 1:length(OGTreat_unique)){
                if(tempID==OGTreat_unique[u]){
                  if(tempID!=OGTreat_unique[k]){
                    temp<-ANOVA_PVALUES[n,]
                    if(length(cup)==0){
                      cup<-temp
                    }
                    else{
                      cup<-rbind(cup,temp)
                    }
                  }
                }
              }
            }
          }

          #alter data to fit with PER,PHA,AMP data
          #make the comparisons easier to read by only listing the one that's not the control
          rows<-cup[,1]
          for (n in 1:length(rows)){
            rows[n] <-paste(".",rows[n],".",sep="")
          }
          spy<-c()
          for(i in 1:length(rows)){
            chonks2 <-rows[i]
            #take out control and separating - so that only the treatment tested against control remains
            test1<-paste("-",OGTreat_unique[k],".",sep="")
            test2<-paste(".",OGTreat_unique[k],"-",sep="")

            if(str_detect(chonks2, test1)){
              chunky<-unlist(strsplit(chonks2,split=test1))
              chunky<-chunky[1]

              #collect the new row labels (treatment names tested against control)
              if(i==1){
                spy<-chunky
              }
              else{
                spy<-rbind(spy,chunky)
              }
            }
            #if it doesn't find one orientation, try the other
            else {
              if(str_detect(chonks2, test2)){
                chunky<-unlist(strsplit(chonks2,split=test2))
                chunky<-chunky[2]

                #collect the new row labels (treatment names tested against control)
                if(i==1){
                  spy<-chunky
                }
                else{
                  spy<-rbind(spy,chunky)
                }
              }
            }
          }
          cup<-cup[,-1]
          #see if marker is still present and remove
          spy<-na.omit(spy)
          for(y in 1:length(spy)){
            if(str_detect(spy[y],".")){
              temp12<-spy[y]
              temp14<-unlist(strsplit(temp12,".",fixed=TRUE))
              for(w in OGTreat_unique){
                if(temp14[1]==w){
                  spy[y]<-temp14[1]
                }
                if(length(temp14)>1){
                  if(temp14[-1]==w){
                    spy[y]<-temp14[-1]
                  }
                }
              }
            }
          }
          spy<-unique(spy)
          #combine labels and data again
          cup<-cbind(spy,cup)

          #need new column names
          hoho<-colnames(my_data)

          colnames(cup)<-hoho
          rownames(cup)<-NULL

          if ((length(hoho)-1)==1){
            P_value<-cbind("P_value","..")
            colnames(P_value)<-hoho
          }

          #make a list of 1's for the control vs self
          one<-c()
          for (o in 1:(length(hoho)-1)){
            one<-cbind(one,"1")
          }
          one<-cbind(as.character(OGTreat_unique[k]),one)
          colnames(one)<-hoho

          #make final file
          if(k==1){
            finale<-rbind(P_value,one,cup)
          }
          else{
            #check if even
            #if even, just add spacer and more p-values
            if(k%%2==0){
              finale<-rbind(finale,P_value,one,cup)
            }
            #check if odd
            #if odd, refresh the labels at top (per,pha,amp) before spacer and more p-values
            if(k%%2!=0){
              finale<-rbind(finale,hoho,P_value,one,cup)
            }
          }
          row.names(finale)<-NULL
          colnames(finale)<-c("Line","ROS_Sum_p_values")
        }

        #create outfile name
        outo <- paste('ANOVA p_values.csv')

        #write outfile
        write.csv(x=finale,outo)

        #move easy read results for ANOVA
        file.move(file=outo, ANOVAFolder, overwrite=TRUE)
        #move DATA for ANOVA
        file.move(file=title2, ANOVAFolder, overwrite=TRUE)
        #remove full results of ANOVA
        file.remove(file=otheroutA)
  }


###################################################
#########End of code###############################
###################################################
