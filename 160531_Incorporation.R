###############################################################################
### Install packages
###############################################################################

if (!require("gridExtra")) {
  install.packages("gridExtra", dependencies = TRUE)
  library(gridExtra)
}

if (!require("ggplot2")) {
  install.packages("ggplot2", dependencies = TRUE)
  library(ggplot2)
}

if (!require("RColorBrewer")) {
  install.packages("RColorBrewer", dependencies = TRUE)
  library(RColorBrewer)
}

if (!require("svDialogs")) {
  install.packages("svDialogs", dependencies = TRUE)
  library(svDialogs)
}

###############################################################################
### Specify SILAC type ###
###############################################################################

inputFolder <- choose.dir(caption = "Select combined/txt/ folder")

projectNameInput <- dlgInput(message = "Name your project:")

SILACType <- select.list(choices = c("LysAndArg", "Leu", "Pro"),
                         graphics = TRUE,
                         title = "Specify labelled amino acids:")

setwd(inputFolder)

projectName <- projectNameInput$res

###############################################################################
### Upload data ###
###############################################################################


peptides <- read.table("peptides.txt",
                       quote = "\"",
                       header = TRUE,
                       sep = "\t",
                       stringsAsFactors = FALSE,
                       comment.char = "")


###############################################################################
### Functions ###
###############################################################################

which.median <- function(x) which.min(abs(x-median(x)))


###############################################################################
### Filter data and transform data ###
###############################################################################
# mark Reverse Hits, Contamintas and OIbS hits as FALSE, others as TRUE
filteredpeptides <- peptides[peptides$Reverse != "+", ]

filteredpeptides <- filteredpeptides[
  filteredpeptides$Potential.contaminant != "+", ]

#filter out peptides with no H/L Ratio
filteredpeptides <- filteredpeptides[
  is.na(filteredpeptides$Ratio.H.L) !=TRUE, ]

#add data for "total" factor
filteredpeptides2 <- filteredpeptides

filteredpeptides2$Peptide <- "Total-Peptides"

if(SILACType == "LysAndArg"){
  #filter for Arginine and Lsy containing peptides
  filteredpeptides$Peptide[filteredpeptides$R.Count > 0 &
                             filteredpeptides$K.Count == 0] <- "Arg-Peptides"
  
  filteredpeptides$Peptide[filteredpeptides$K.Count > 0 &
                             filteredpeptides$R.Count == 0] <- "Lys-Peptides"
  
  filteredpeptides <- filteredpeptides[!is.na(filteredpeptides$Peptide), ]
  
  #merge total and AS specific data
  filteredpeptides <- rbind(filteredpeptides, filteredpeptides2)
  
  # Calculate incorporation rate per peptide
  
  filteredpeptides$IncorporationRate <- 1 - 1 / (filteredpeptides$Ratio.H.L+1)
  
  DenK <- density(1-1/(filteredpeptides$Ratio.H.L[
    filteredpeptides$Peptide == "Lys-Peptides"]+1))
  
  DenR <- density(1-1/(filteredpeptides$Ratio.H.L[
    filteredpeptides$Peptide == "Arg-Peptides"]+1))
  
  DenAll <- density(1-1/(filteredpeptides$Ratio.H.L[
    filteredpeptides$Peptide == "Total-Peptides"]+1))
  
  IncK <- round(x = DenK$x[which.max(DenK$y)], digits = 4)
  
  IncR <- round(x = DenR$x[which.max(DenR$y)], digits = 4)
  
  IncAll <- round(x = DenAll$x[which.max(DenAll$y)], digits = 4)
  
  MedianK <- round(x = median(filteredpeptides$IncorporationRate[
    filteredpeptides$Peptide == "Lys-Peptides"]), digits = 4)
  
  MedianR <- round(x = median(filteredpeptides$IncorporationRate[
    filteredpeptides$Peptide == "Arg-Peptides"]), digits = 4)
  
  MedianAll <- round(x = median(filteredpeptides$IncorporationRate[
    filteredpeptides$Peptide == "Total-Peptides"]), digits = 4)
  
  Peptides <- c("Lys-Peptides", "Arg-Peptides", "Total-Peptides")
  
  Maximum.Density <- c(IncK, IncR, IncAll)
  
  Number <- c(length(
    filteredpeptides$Peptide[
      filteredpeptides$Peptide == "Lys-Peptides"]),
              length(
                filteredpeptides$Peptide[
                  filteredpeptides$Peptide == "Arg-Peptides"]),
              length(
                filteredpeptides$Peptide[
                  filteredpeptides$Peptide == "Total-Peptides"]))
  
  Median.Density <- c(MedianK, MedianR, MedianAll)
  
  myTable <- data.frame(Peptides, Number, Maximum.Density)
  
  maxY <- max(max(DenK$y), max(DenR$y), max(DenAll$y))
  
} else {
  
  if(SILACType == "Leu"){
    
    #filter for Leu containing peptides
    filteredpeptides$Peptide[filteredpeptides$L.Count > 0] <- "Leu-Peptides"
    
    filteredpeptides <- filteredpeptides[!is.na(filteredpeptides$Peptide), ]
    
    # Calculate incorporation rate per peptide
    
    filteredpeptides$IncorporationRate <- 1 - 1 / (
      filteredpeptides$Ratio.H.L+1)
    
    DenL <- density(1-1/filteredpeptides$Ratio.H.L)
    
    IncL <- round(x = DenL$x[which.max(DenL$y)], digits = 4)
    
    MedianL <- round(x = median(filteredpeptides$IncorporationRate[
      filteredpeptides$Peptide == "Leu-Peptides"]), digits = 4)
    
    Peptides <- "Leu-Peptides"
    
    Maximum.Density <- IncL
    
    Number <- length(filteredpeptides$Peptide)
    
    myTable <- data.frame(Peptides, Number, Maximum.Density)
    
  } else {
    
    if(SILACType == "Pro"){
      
      #filter for Pro containing peptides
      filteredpeptides$Peptide[filteredpeptides$P.Count > 0] <- "Pro-Peptides"
      
      filteredpeptides <- filteredpeptides[!is.na(filteredpeptides$Peptide), ]
      
      # Calculate incorporation rate per peptide
      
      filteredpeptides$IncorporationRate <- 1 - 1 / (
        filteredpeptides$Ratio.H.L+1)
      
      DenP <- density(1-1/filteredpeptides$Ratio.H.L)
      
      IncP <- round(x = DenP$x[which.max(DenP$y)], digits = 4)
      
      MedianP <- round(x = median(filteredpeptides$IncorporationRate[
        filteredpeptides$Peptide == "Pro-Peptides"]), digits = 4)
      
      Peptides <- "Pro-Peptides"
      
      Maximum.Density <- IncP
      
      Number <- length(filteredpeptides$Peptide)
      
      myTable <- data.frame(Peptides, Number, Maximum.Density)
      
    }
    
  }
  
}


#Create incorporation plot and export as pdf

IncPlot <- ggplot(data = filteredpeptides, aes(x = IncorporationRate,
                                               fill = Peptide,
                                               colour = Peptide)) +
  
  geom_density(alpha = 0.3) +
  
  scale_fill_manual(values = c(brewer.pal(8,"Dark2"))) +
  
  scale_colour_manual(values = c(brewer.pal(8,"Dark2"))) +
  
  scale_x_continuous(limits = c(0, 1)) +
  
  geom_vline(xintercept = 0.95, linetype = "dashed", colour = "red") +
  
  guides(alpha = FALSE) +
  
  annotation_custom(tableGrob(myTable, rows = NULL),
                    xmin = 0.1,
                    xmax = 0.2) +
  
  geom_text(aes(0.945, 0, label = "95 % cuttoff"),
            colour ="red",
            hjust = 1,
            vjust = 0)+
  
  ggtitle(paste0("Silac Incorporation rate project: ", projectName))

pdf(paste0(projectName, "_", "SILAC_incorp_rate.pdf"),
    width = 20.6,
    height = 11.6)

IncPlot

dev.off()
