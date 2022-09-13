
#Plotting the coverage of one chromosome (or mitogenome in my case) in R #####
#code: 
library(tidyverse)

samtools depth ${prefix}_${base}_sorted.bam > deduped_${prefix}_${base}.coverage

#I know that ascaris worked best, however, I will need to adapt this for multiple chromosomes????
awk '$1 == "NC_016198_Ascaris_lumbricoides_mitochondrion_complete_genome" {print $0}' deduped_TKU102.coverage > ascaris_${prefix}_${base}.coverage


#calculate the depth per sample 
samtools depth ijp_human_mito_TKU102_trimmed_sorted.bam > deduped_TKU102.coverage 

#isolate the results only for ascaris 
awk '$1 == "NC_016198_Ascaris_lumbricoides_mitochondrion_complete_genome" {print $0}' deduped_TKU102.coverage > ascarisTKU102.coverage


#take the file ".coverage" into R and plot it 
library(tidyverse)
---------------------------
#TKU102 ----
TKU102_ascaris <- read.table("/Users/marinapapaiakovou/Documents/00.Cambridge_PhD/02.Science/02.Genome_skimming/2019_Metagenomics_NHM/IJP_genome_skimming/Data/ascarisTKU102.coverage",
                      header=FALSE, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)

library(reshape) # loads the library to rename the column names

TKU102_ascaris <- plyr::rename(TKU102_ascaris,c(V1="reference", V2="locus", V3="depth")) # renames the header

#Now, plot the coverage by depth:
  
#plot(TKU102_ascaris$locus, TKU102_ascaris$depth)

#to get the wider/cleaner view of the plot use

library(lattice, pos=10) 
plot1 <- xyplot(depth ~ locus, type="p", pch=20, auto.key=list(border=TRUE), par.settings=simpleTheme(col="black",pch=20), 
                                scales=list(x=list(relation='same'), y=list(relation='same')), 
                                data=TKU102_ascaris, main="depth by locus - Ascaris lumbricoides (TKU102)")

---------------------------------------------------------
  
#NDK63-----  
#trying to plot another BAM file and then put them in a single plot 
NDK63_ascaris <- read.table("/Users/marinapapaiakovou/Documents/00.Cambridge_PhD/02.Science/02.Genome_skimming/2019_Metagenomics_NHM/IJP_genome_skimming/Data/ascarisNDK63.coverage",
                             header=FALSE, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)

library(reshape) # loads the library to rename the column names

NDK63_ascaris <- plyr::rename(NDK63_ascaris,c(V1="reference", V2="locus", V3="depth")) # renames the header

#Now, plot the coverage by depth:

plot(all$locus, all$depth)

#to get the wider/cleaner view of the plot use

library(lattice, pos=10) 
plot2 <- xyplot(depth ~ locus, type="p", pch=20, auto.key=list(border=TRUE), par.settings=simpleTheme(col="black",pch=20), 
       scales=list(x=list(relation='same'), y=list(relation='same')), 
       data=NDK63_ascaris, main="depth by locus - Ascaris lumbricoides (NDK63)")

---------------------------
#used to have ET018
-----------------------------------------

#ET103-----  
#trying to plot another BAM file and then put them in a single plot 
ET103_ascaris <- read.table("/Users/marinapapaiakovou/Documents/00.Cambridge_PhD/02.Science/02.Genome_skimming/2019_Metagenomics_NHM/IJP_genome_skimming/Data/ascarisET103.coverage",
                            header=FALSE, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)

library(reshape) # loads the library to rename the column names

ET103_ascaris <- plyr::rename(ET103_ascaris,c(V1="reference", V2="locus", V3="depth")) # renames the header

#Now, plot the coverage by depth:

plot(all$locus, all$depth)

#to get the wider/cleaner view of the plot use

library(lattice, pos=10) 
plot4 <- xyplot(depth ~ locus, type="p", pch=20, auto.key=list(border=TRUE), par.settings=simpleTheme(col="black",pch=20), 
                scales=list(x=list(relation='same'), y=list(relation='same')), 
                data=ET103_ascaris, main="depth by locus - Ascaris lumbricoides (ET103)")



-----------------------------------------
#NDK113 -----
NDK113_ascaris <- read.table("/Users/marinapapaiakovou/Documents/00.Cambridge_PhD/02.Science/02.Genome_skimming/2019_Metagenomics_NHM/IJP_genome_skimming/Data/ascarisNDK113.coverage",
                            header=FALSE, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)

library(reshape) # loads the library to rename the column names

NDK113_ascaris <- plyr::rename(NDK113_ascaris,c(V1="reference", V2="locus", V3="depth")) # renames the header

#Now, plot the coverage by depth:

plot(all$locus, all$depth)

#to get the wider/cleaner view of the plot use

library(lattice, pos=10) 
plot5 <- xyplot(depth ~ locus, type="p", pch=20, auto.key=list(border=TRUE), par.settings=simpleTheme(col="black",pch=20), 
                scales=list(x=list(relation='same'), y=list(relation='same')), 
                data=NDK113_ascaris, main="depth by locus - Ascaris lumbricoides (NDK113)")

  
-----------------------------------------  

#NDK92 ----  

NDK92_ascaris <- read.table("/Users/marinapapaiakovou/Documents/00.Cambridge_PhD/02.Science/02.Genome_skimming/2019_Metagenomics_NHM/IJP_genome_skimming/Data/ascarisNDK92.coverage",
                               header=FALSE, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)

library(reshape) # loads the library to rename the column names

NDK92_ascaris <- plyr::rename(NDK92_ascaris,c(V1="reference", V2="locus", V3="depth")) # renames the header

#Now, plot the coverage by depth:

plot(all$locus, all$depth)

#to get the wider/cleaner view of the plot use

library(lattice, pos=10) 
plot6<- xyplot(depth ~ locus, type="p", pch=20, auto.key=list(border=TRUE), par.settings=simpleTheme(col="black",pch=20), 
                scales=list(x=list(relation='same'), y=list(relation='same')), 
                data=NDK92_ascaris, main="depth by locus - Ascaris lumbricoides (NDK92)")  

-----------------------------------------  

#TKU23----
TKU23_ascaris <- read.table("/Users/marinapapaiakovou/Documents/00.Cambridge_PhD/02.Science/02.Genome_skimming/2019_Metagenomics_NHM/IJP_genome_skimming/Data/ascarisTKU23.coverage",
                              header=FALSE, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)

library(reshape) # loads the library to rename the column names

TKU23_ascaris  <- plyr::rename(TKU23_ascaris ,c(V1="reference", V2="locus", V3="depth")) # renames the header

#Now, plot the coverage by depth:

plot(all$locus, all$depth)

#to get the wider/cleaner view of the plot use

library(lattice, pos=10) 
plot7<- xyplot(depth ~ locus, type="p", pch=20, auto.key=list(border=TRUE), par.settings=simpleTheme(col="black",pch=20), 
               scales=list(x=list(relation='same'), y=list(relation='same')), 
               data=TKU23_ascaris , main="depth by locus - Ascaris lumbricoides (TKU23)")  

----------------------------------------    
#TKU25 ----

TKU25_ascaris <- read.table("/Users/marinapapaiakovou/Documents/00.Cambridge_PhD/02.Science/02.Genome_skimming/2019_Metagenomics_NHM/IJP_genome_skimming/Data/ascarisTKU25.coverage",
                            header=FALSE, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)

library(reshape) # loads the library to rename the column names

TKU25_ascaris <- plyr::rename(TKU25_ascaris,c(V1="reference", V2="locus", V3="depth")) # renames the header

#Now, plot the coverage by depth:

plot(all$locus, all$depth)

#to get the wider/cleaner view of the plot use

library(lattice, pos=10) 
plot8<- xyplot(depth ~ locus, type="p", pch=20, auto.key=list(border=TRUE), par.settings=simpleTheme(col="black",pch=20), 
               scales=list(x=list(relation='same'), y=list(relation='same')), 
               data=TKU25_ascaris, main="depth by locus - Ascaris lumbricoides (TKU25)")  
-----------------------------------------        
#plot the two samples together 
install.packages("gridExtra")
library(gridExtra)
grid.arrange(plot1,plot2, plot4, nrow=3)

grid.arrange(plot5, plot6, plot7,plot8, nrow=4)
