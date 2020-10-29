############################################################################
# Barn Swallow (China) mate-pair variable correlations
############################################################################

### Goal: examine preliminary relationships between isotope, phenotype, and
### ancestry variables for Barn Swallow mate-pairs from the China rustica x
### gutturalis hybrid zone.

### Clear environment-------------------------------------------------------

rm(list = ls())

### Load dependencies-------------------------------------------------------

install.packages("rstatix")
install.packages("ggpubr")
install.packages("tibble")

library(Rmisc)
library(rstatix)
library(ggpubr)

### Read in data------------------------------------------------------------

## See associated README file for the Barn Swallow assortativity project for
## details on the code used to relate non-redundant pairs. The convention will
## be that the variable for the second mate will be the same as the first, with
## ".1" added to the string (e.g., mate 1 = d2h; mate 2 = d2h.1).

mates <- read.table('/Volumes/GoogleDrive/My Drive/Drew_and_Becca_shared_folder/projects/hirundo_assortativity/data/PAIRS_phenotype+ancestry+isotope_assignment+geo_data_08.10.20.FIX.noGPS.txt',header=T)

### Set up color palette to color points by locations------------------------

# By location: Jiuquan, Lanzhou, Zhangye
loc.pal <- palette(c('red','turquoise3','orange'))

### Plot scatterplots and test correlations between variables for mates-----

## Isotopes

par(mfrow=c(1,2))
plot(mates$d2h,mates$d2h.1,pch=20,col=mates$location)
plot(mates$d13C_VPDB,mates$d13C_VPDB.1,pch=20,col=mates$location)
legend("topleft",legend=unique(mates$location),col=1:length(mates$location),pch=20)


# Are isotope variables correlated for mates?
cor.test(mates$d2h,mates$d2h.1,method="spearman")
# Yes. ------------------------------------------------------------------------------------------***
cor.test(mates$d13C_VPDB,mates$d13C_VPDB.1,method="spearman")
# No.

## Ancestry

par(mfrow=c(1,2))
plot(mates$proportion.rustica.ancestry,mates$proportion.rustica.ancestry.1,pch=20,col=mates$location)
plot(mates$proportion.gutturalis.ancestry,mates$proportion.gutturalis.ancestry.1,pch=20,col=mates$location)
legend("topleft",legend=unique(mates$location),col=1:length(mates$location),pch=20)

# Are ancestry coefficients correlated?
cor.test(mates$proportion.rustica.ancestry,mates$proportion.rustica.ancestry.1,method="spearman")
# Yes. --------------------------------------------------------------------------------------------------------------------***
cor.test(mates$proportion.gutturalis.ancestry,mates$proportion.gutturalis.ancestry.1,method="spearman")
# Yes. Just the inverse ancestry coefficients from above.

## Wing/tail streamer lengths

par(mfrow=c(1,3))
plot(mates$mean_rwl,mates$mean_rwl.1,pch=20,col=mates$location)
plot(mates$mean_rts,mates$mean_rts.1,pch=20,xlim=c(70,130),ylim=c(70,130),col=mates$location)
plot(mates$mean_lts,mates$mean_lts.1,pch=20,xlim=c(70,130),ylim=c(70,130),col=mates$location)
legend("topleft",legend=unique(mates$location),col=1:length(mates$location),pch=20)

# Are wing/tail streamer lengths correlated?
cor.test(mates$mean_rwl,mates$mean_rwl.1,method="spearman")
# No.
cor.test(mates$mean_rts,mates$mean_rts.1,method="spearman")
# No. 
cor.test(mates$mean_lts,mates$mean_lts.1,method="spearman")
# Yes, positive. ------------------------------------------------------------------------------------------***

## Belly feathers

par(mfrow=c(1,4))
plot(mates$belly.total.bright,mates$belly.total.bright.1,pch=20,col=mates$location)
plot(mates$belly.avg.bright,mates$belly.avg.bright.1,pch=20,col=mates$location)
plot(mates$belly.hue,mates$belly.hue.1,pch=20,col=mates$location)
plot(mates$belly.chrom,mates$belly.chrom.1,pch=20,col=mates$location)
legend("topleft",legend=unique(mates$location),col=1:length(mates$location),pch=20)

# Are belly feather variables correlated?
cor.test(mates$belly.total.bright,mates$belly.total.bright.1,method="spearman")
# No.
cor.test(mates$belly.avg.bright,mates$belly.avg.bright.1,method="spearman")
# No.
cor.test(mates$belly.hue,mates$belly.hue.1,method="spearman")
# Yes. ------------------------------------------------------------------------------------------***
cor.test(mates$belly.chrom,mates$belly.chrom.1,method="spearman")
# No.

## Throat feathers

par(mfrow=c(1,4))
plot(mates$throat.total.bright,mates$throat.total.bright.1,pch=20,col=mates$location)
plot(mates$throat.avg.bright,mates$throat.avg.bright.1,pch=20,col=mates$location)
plot(mates$throat.hue,mates$throat.hue.1,pch=20,col=mates$location)
plot(mates$throat.chrom,mates$throat.chrom.1,pch=20,col=mates$location)
legend("topleft",legend=unique(mates$location),col=1:length(mates$location),pch=20)

# Are throat feather variables correlated?
cor.test(mates$throat.total.bright,mates$throat.total.bright.1,method="spearman")
# YES.------------------------------------------------------------------------------------------***
cor.test(mates$throat.avg.bright,mates$throat.avg.bright.1,method="spearman")
# YES. ------------------------------------------------------------------------------------------***
cor.test(mates$throat.hue,mates$throat.hue.1,method="spearman")
# No.
cor.test(mates$throat.chrom,mates$throat.chrom.1,method="spearman")
# No.

## Breast feathers

par(mfrow=c(1,4))
plot(mates$breast.total.bright,mates$breast.total.bright.1,pch=20,col=mates$location)
plot(mates$breast.avg.bright,mates$breast.avg.bright.1,pch=20,col=mates$location)
plot(mates$breast.hue,mates$breast.hue.1,pch=20,col=mates$location)
plot(mates$breast.chrom,mates$breast.chrom.1,pch=20,col=mates$location)
legend("topleft",legend=unique(mates$location),col=1:length(mates$location),pch=20)

# Are breast feather variables correlated?
cor.test(mates$breast.total.bright,mates$breast.total.bright.1,method="spearman")
# No.
cor.test(mates$breast.avg.bright,mates$breast.avg.bright.1,method="spearman")
# No.
cor.test(mates$breast.hue,mates$breast.hue.1,method="spearman")
# YES.------------------------------------------------------------------------------------------***
cor.test(mates$breast.chrom,mates$breast.chrom.1,method="spearman")
# No.

## Vent feathers

par(mfrow=c(1,4))
plot(mates$vent.total.bright,mates$vent.total.bright.1,pch=20,col=mates$location)
plot(mates$vent.avg.bright,mates$vent.avg.bright.1,pch=20,col=mates$location)
plot(mates$vent.hue,mates$vent.hue.1,pch=20,col=mates$location)
plot(mates$vent.chrom,mates$vent.chrom.1,pch=20,col=mates$location)
legend("topleft",legend=unique(mates$location),col=1:length(mates$location),pch=20)

# Are vent feather variables correlated?
cor.test(mates$vent.total.bright,mates$vent.total.bright.1,method="spearman")
# No.
cor.test(mates$vent.avg.bright,mates$vent.avg.bright.1,method="spearman")
# No.
cor.test(mates$vent.hue,mates$vent.hue.1,method="spearman")
# Yes. ------------------------------------------------------------------------------------------***
cor.test(mates$vent.chrom,mates$vent.chrom.1,method="spearman")
# No.

## Bird mass

par(mfrow=c(1,1))
plot(mates$bird_mass,mates$bird_mass.1,pch=20,col=mates$location)
legend("topleft",legend=unique(mates$location),col=1:length(mates$location),pch=20)

# Is bird mass correlated?
cor.test(mates$bird_mass,mates$bird_mass.1,method="spearman")
# Yes. ------------------------------------------------------------------------------------------***

plot(mates$bird_mass,mates$proportion.rustica.ancestry,pch=20,col=mates$location)
cor.test(mates$bird_mass,mates$proportion.rustica.ancestry,method="spearman")

### -------------------------------------------------------------------------
### Analysis of overall relationships
### -------------------------------------------------------------------------

### Compute correlation matrix for all variables

mates.prune <- mates[, c(16,17,21,22,23,25,29,33,37,51,70,87,88,92,93,94,96,100,104,108,122,141)]
head(mates.prune, 5)

# Set correlation matrix
mates.cor_mat <- mates.prune %>% cor_mat()
# Get p-values
mates.cor_mat %>% cor_get_pval()

write.table(mates.cor_mat, "/Volumes/GoogleDrive/My Drive/github/assortativity/data/correlation_table.txt",quote=FALSE,row.names=FALSE,sep='\t')
write.table(mates.cor_mat %>% cor_get_pval(), "/Volumes/GoogleDrive/My Drive/github/assortativity/data/correlation_table_p-val.txt",quote=FALSE,row.names=FALSE,sep='\t')

# Specify color palette and plot correlation matrix
cor.palette <- get_palette("PuOr",200)
china.cor_mat %>%
  pull_lower_triangle() %>%
  cor_plot(label = FALSE,insignificant='cross',method='circle',palette = cor.palette)
# Output long-format matrix
write.table(china.cor_mat %>% cor_gather(),"correlation_table_long.txt",quote=FALSE,row.names=FALSE,sep='\t')


### ------------------------------------------------------------------------
### Examine relationships between male and female variables-----------------
### ------------------------------------------------------------------------

### Here, the goal is to look at preliminary Spearman's correlations between
### paired male and female variables, looking in both directions...

# FEMALE variables against MALE ancestry

par(mfrow=c(4,3))

plot(mates$proportion.rustica.ancestry,mates$d2h.1,pch=20,col=mates$location)
plot(mates$proportion.rustica.ancestry,mates$d13C_VPDB.1,pch=20,col=mates$location)
plot(mates$proportion.rustica.ancestry,mates$mean_rwl.1,pch=20,col=mates$location)
plot(mates$proportion.rustica.ancestry,mates$mean_rts.1,pch=20,col=mates$location)
plot(mates$proportion.rustica.ancestry,mates$mean_lts.1,pch=20,col=mates$location)
plot(mates$proportion.rustica.ancestry,mates$belly.avg.bright.1,pch=20,col=mates$location)
plot(mates$proportion.rustica.ancestry,mates$breast.avg.bright.1,pch=20,col=mates$location)
plot(mates$proportion.rustica.ancestry,mates$throat.avg.bright.1,pch=20,col=mates$location)
plot(mates$proportion.rustica.ancestry,mates$vent.avg.bright.1,pch=20,col=mates$location)
plot(mates$proportion.rustica.ancestry,mates$breast.avg.bright.1,pch=20,col=mates$location)
plot(mates$proportion.rustica.ancestry,mates$bird_mass.1,pch=20,col=mates$location)

cor.test(mates$proportion.rustica.ancestry,mates$d2h.1,method="spearman")                               ###
cor.test(mates$proportion.rustica.ancestry,mates$d13C_VPDB.1,method="spearman")                         ###
cor.test(mates$proportion.rustica.ancestry,mates$mean_rwl.1,method="spearman")                          ###
cor.test(mates$proportion.rustica.ancestry,mates$mean_rts.1,method="spearman")
cor.test(mates$proportion.rustica.ancestry,mates$mean_lts.1,method="spearman")
cor.test(mates$proportion.rustica.ancestry,mates$belly.avg.bright.1,method="spearman")
cor.test(mates$proportion.rustica.ancestry,mates$breast.avg.bright.1,method="spearman")
cor.test(mates$proportion.rustica.ancestry,mates$throat.avg.bright.1,method="spearman")
cor.test(mates$proportion.rustica.ancestry,mates$vent.avg.bright.1,method="spearman")                   ###
cor.test(mates$proportion.rustica.ancestry,mates$bird_mass.1,method="spearman")

# MALE variables against FEMALE ancestry

par(mfrow=c(4,3))

plot(mates$proportion.rustica.ancestry.1,mates$d2h,pch=20,col=mates$location)
plot(mates$proportion.rustica.ancestry.1,mates$d13C_VPDB,pch=20,col=mates$location)
plot(mates$proportion.rustica.ancestry.1,mates$mean_rwl,pch=20,col=mates$location)
plot(mates$proportion.rustica.ancestry.1,mates$mean_rts,pch=20,col=mates$location)
plot(mates$proportion.rustica.ancestry.1,mates$mean_lts,pch=20,col=mates$location)
plot(mates$proportion.rustica.ancestry.1,mates$belly.avg.bright,pch=20,col=mates$location)
plot(mates$proportion.rustica.ancestry.1,mates$breast.avg.bright,pch=20,col=mates$location)
plot(mates$proportion.rustica.ancestry.1,mates$throat.avg.bright,pch=20,col=mates$location)
plot(mates$proportion.rustica.ancestry.1,mates$vent.avg.bright,pch=20,col=mates$location)
plot(mates$proportion.rustica.ancestry.1,mates$breast.avg.bright,pch=20,col=mates$location)
plot(mates$proportion.rustica.ancestry.1,mates$bird_mass,pch=20,col=mates$location)

cor.test(mates$proportion.rustica.ancestry.1,mates$d2h,method="spearman")                               
cor.test(mates$proportion.rustica.ancestry.1,mates$d13C_VPDB,method="spearman")                         
cor.test(mates$proportion.rustica.ancestry.1,mates$mean_rwl,method="spearman")                          ###
cor.test(mates$proportion.rustica.ancestry.1,mates$mean_rts,method="spearman")                          ###
cor.test(mates$proportion.rustica.ancestry.1,mates$mean_lts,method="spearman")                          ###
cor.test(mates$proportion.rustica.ancestry.1,mates$belly.avg.bright,method="spearman")
cor.test(mates$proportion.rustica.ancestry.1,mates$breast.avg.bright,method="spearman")
cor.test(mates$proportion.rustica.ancestry.1,mates$throat.avg.bright,method="spearman")
cor.test(mates$proportion.rustica.ancestry.1,mates$vent.avg.bright,method="spearman")                   
cor.test(mates$proportion.rustica.ancestry.1,mates$bird_mass,method="spearman")


# FEMALE variables against MALE hydrogen isotope

par(mfrow=c(4,3))

plot(mates$d2h,mates$proportion.rustica.ancestry.1,pch=20,col=mates$location)
plot(mates$d2h,mates$d13C_VPDB.1,pch=20,col=mates$location)
plot(mates$d2h,mates$mean_rwl.1,pch=20,col=mates$location)
plot(mates$d2h,mates$mean_rts.1,pch=20,col=mates$location)
plot(mates$d2h,mates$mean_lts.1,pch=20,col=mates$location)
plot(mates$d2h,mates$belly.avg.bright.1,pch=20,col=mates$location)
plot(mates$d2h,mates$breast.avg.bright.1,pch=20,col=mates$location)
plot(mates$d2h,mates$throat.avg.bright.1,pch=20,col=mates$location)
plot(mates$d2h,mates$vent.avg.bright.1,pch=20,col=mates$location)
plot(mates$d2h,mates$bird_mass.1,pch=20,col=mates$location)

cor.test(mates$d2h,mates$proportion.rustica.ancestry.1,method="spearman")
cor.test(mates$d2h,mates$d13C_VPDB.1,method="spearman")                         ###
cor.test(mates$d2h,mates$mean_rwl.1,method="spearman")                          
cor.test(mates$d2h,mates$mean_rts.1,method="spearman")
cor.test(mates$d2h,mates$mean_lts.1,method="spearman")
cor.test(mates$d2h,mates$belly.avg.bright.1,method="spearman")                  ###
cor.test(mates$d2h,mates$breast.avg.bright.1,method="spearman")                 ###
cor.test(mates$d2h,mates$throat.avg.bright.1,method="spearman")
cor.test(mates$d2h,mates$vent.avg.bright.1,method="spearman")                   ###
cor.test(mates$d2h,mates$bird_mass.1,method="spearman")

# MALE variables against FEMALE hydrogen isotope

par(mfrow=c(4,3))

plot(mates$d2h.1,mates$proportion.rustica.ancestry,pch=20,col=mates$location)
plot(mates$d2h.1,mates$d13C_VPDB,pch=20,col=mates$location)
plot(mates$d2h.1,mates$mean_rwl,pch=20,col=mates$location)
plot(mates$d2h.1,mates$mean_rts,pch=20,col=mates$location)
plot(mates$d2h.1,mates$mean_lts,pch=20,col=mates$location)
plot(mates$d2h.1,mates$belly.avg.bright,pch=20,col=mates$location)
plot(mates$d2h.1,mates$breast.avg.bright,pch=20,col=mates$location)
plot(mates$d2h.1,mates$throat.avg.bright,pch=20,col=mates$location)
plot(mates$d2h.1,mates$vent.avg.bright,pch=20,col=mates$location)
plot(mates$d2h.1,mates$bird_mass,pch=20,col=mates$location)

cor.test(mates$d2h.1,mates$proportion.rustica.ancestry,method="spearman")       ###
cor.test(mates$d2h.1,mates$d13C_VPDB,method="spearman")                        
cor.test(mates$d2h.1,mates$mean_rwl,method="spearman")                         
cor.test(mates$d2h.1,mates$mean_rts,method="spearman")
cor.test(mates$d2h.1,mates$mean_lts,method="spearman")
cor.test(mates$d2h.1,mates$belly.avg.bright,method="spearman")                 
cor.test(mates$d2h.1,mates$breast.avg.bright,method="spearman")                
cor.test(mates$d2h.1,mates$throat.avg.bright,method="spearman")
cor.test(mates$d2h.1,mates$vent.avg.bright,method="spearman")                  
cor.test(mates$d2h.1,mates$bird_mass,method="spearman")


# FEMALE variables against MALE carbon isotope

par(mfrow=c(4,3))

plot(mates$d13C_VPDB,mates$proportion.rustica.ancestry.1,pch=20,col=mates$location)
plot(mates$d13C_VPDB,mates$d2h.1,pch=20,col=mates$location)
plot(mates$d13C_VPDB,mates$mean_rwl.1,pch=20,col=mates$location)
plot(mates$d13C_VPDB,mates$mean_rts.1,pch=20,col=mates$location)
plot(mates$d13C_VPDB,mates$mean_lts.1,pch=20,col=mates$location)
plot(mates$d13C_VPDB,mates$belly.avg.bright.1,pch=20,col=mates$location)
plot(mates$d13C_VPDB,mates$breast.avg.bright.1,pch=20,col=mates$location)
plot(mates$d13C_VPDB,mates$throat.avg.bright.1,pch=20,col=mates$location)
plot(mates$d13C_VPDB,mates$vent.avg.bright.1,pch=20,col=mates$location)
plot(mates$d13C_VPDB,mates$bird_mass.1,pch=20,col=mates$location)

cor.test(mates$d13C_VPDB,mates$proportion.rustica.ancestry.1,method="spearman")
cor.test(mates$d13C_VPDB,mates$d2h.1,method="spearman")                         
cor.test(mates$d13C_VPDB,mates$mean_rwl.1,method="spearman")                          
cor.test(mates$d13C_VPDB,mates$mean_rts.1,method="spearman")
cor.test(mates$d13C_VPDB,mates$mean_lts.1,method="spearman")
cor.test(mates$d13C_VPDB,mates$belly.avg.bright.1,method="spearman")                  
cor.test(mates$d13C_VPDB,mates$breast.avg.bright.1,method="spearman")                 ###
cor.test(mates$d13C_VPDB,mates$throat.avg.bright.1,method="spearman")                 ###
cor.test(mates$d13C_VPDB,mates$vent.avg.bright.1,method="spearman")                   
cor.test(mates$d13C_VPDB,mates$bird_mass.1,method="spearman")

# MALE variables against FEMALE carbon isotope

par(mfrow=c(4,3))

plot(mates$d13C_VPDB.1,mates$proportion.rustica.ancestry,pch=20,col=mates$location)
plot(mates$d13C_VPDB.1,mates$d2h,pch=20,col=mates$location)
plot(mates$d13C_VPDB.1,mates$mean_rwl,pch=20,col=mates$location)
plot(mates$d13C_VPDB.1,mates$mean_rts,pch=20,col=mates$location)
plot(mates$d13C_VPDB.1,mates$mean_lts,pch=20,col=mates$location)
plot(mates$d13C_VPDB.1,mates$belly.avg.bright,pch=20,col=mates$location)
plot(mates$d13C_VPDB.1,mates$breast.avg.bright,pch=20,col=mates$location)
plot(mates$d13C_VPDB.1,mates$throat.avg.bright,pch=20,col=mates$location)
plot(mates$d13C_VPDB.1,mates$vent.avg.bright,pch=20,col=mates$location)
plot(mates$d13C_VPDB.1,mates$bird_mass,pch=20,col=mates$location)

cor.test(mates$d13C_VPDB.1,mates$proportion.rustica.ancestry,method="spearman")       ###
cor.test(mates$d13C_VPDB.1,mates$d2h,method="spearman")                               ###                 
cor.test(mates$d13C_VPDB.1,mates$mean_rwl,method="spearman")                          ###            
cor.test(mates$d13C_VPDB.1,mates$mean_rts,method="spearman")
cor.test(mates$d13C_VPDB.1,mates$mean_lts,method="spearman")
cor.test(mates$d13C_VPDB.1,mates$belly.avg.bright,method="spearman")                 
cor.test(mates$d13C_VPDB.1,mates$breast.avg.bright,method="spearman")                
cor.test(mates$d13C_VPDB.1,mates$throat.avg.bright,method="spearman")
cor.test(mates$d13C_VPDB.1,mates$vent.avg.bright,method="spearman")                  
cor.test(mates$d13C_VPDB.1,mates$bird_mass,method="spearman")


# FEMALE variables against MALE right wing length

par(mfrow=c(4,3))

plot(mates$mean_rwl,mates$proportion.rustica.ancestry.1,pch=20,col=mates$location)
plot(mates$mean_rwl,mates$d2h.1,pch=20,col=mates$location)
plot(mates$mean_rwl,mates$d13C_VPDB.1,pch=20,col=mates$location)
plot(mates$mean_rwl,mates$mean_rts.1,pch=20,col=mates$location)
plot(mates$mean_rwl,mates$mean_lts.1,pch=20,col=mates$location)
plot(mates$mean_rwl,mates$belly.avg.bright.1,pch=20,col=mates$location)
plot(mates$mean_rwl,mates$breast.avg.bright.1,pch=20,col=mates$location)
plot(mates$mean_rwl,mates$throat.avg.bright.1,pch=20,col=mates$location)
plot(mates$mean_rwl,mates$vent.avg.bright.1,pch=20,col=mates$location)
plot(mates$mean_rwl,mates$bird_mass.1,pch=20,col=mates$location)

cor.test(mates$mean_rwl,mates$proportion.rustica.ancestry.1,method="spearman")       ###
cor.test(mates$mean_rwl,mates$d2h.1,method="spearman")                         
cor.test(mates$mean_rwl,mates$d13C_VPDB.1,method="spearman")                         ###                   
cor.test(mates$mean_rwl,mates$mean_rts.1,method="spearman")
cor.test(mates$mean_rwl,mates$mean_lts.1,method="spearman")
cor.test(mates$mean_rwl,mates$belly.avg.bright.1,method="spearman")                  
cor.test(mates$mean_rwl,mates$breast.avg.bright.1,method="spearman")                 ###
cor.test(mates$mean_rwl,mates$throat.avg.bright.1,method="spearman")                 
cor.test(mates$mean_rwl,mates$vent.avg.bright.1,method="spearman")                   ###        
cor.test(mates$mean_rwl,mates$bird_mass.1,method="spearman")                         ###

# MALE variables against FEMALE right wing length

par(mfrow=c(4,3))

plot(mates$mean_rwl.1,mates$proportion.rustica.ancestry,pch=20,col=mates$location)
plot(mates$mean_rwl.1,mates$d2h,pch=20,col=mates$location)
plot(mates$mean_rwl.1,mates$d13C_VPDB,pch=20,col=mates$location)
plot(mates$mean_rwl.1,mates$mean_rts,pch=20,col=mates$location)
plot(mates$mean_rwl.1,mates$mean_lts,pch=20,col=mates$location)
plot(mates$mean_rwl.1,mates$belly.avg.bright,pch=20,col=mates$location)
plot(mates$mean_rwl.1,mates$breast.avg.bright,pch=20,col=mates$location)
plot(mates$mean_rwl.1,mates$throat.avg.bright,pch=20,col=mates$location)
plot(mates$mean_rwl.1,mates$vent.avg.bright,pch=20,col=mates$location)
plot(mates$mean_rwl.1,mates$bird_mass,pch=20,col=mates$location)

cor.test(mates$mean_rwl.1,mates$proportion.rustica.ancestry,method="spearman")       ###
cor.test(mates$mean_rwl.1,mates$d2h,method="spearman")                                             
cor.test(mates$mean_rwl.1,mates$d13C_VPDB,method="spearman")                                    
cor.test(mates$mean_rwl.1,mates$mean_rts,method="spearman")                          ###
cor.test(mates$mean_rwl.1,mates$mean_lts,method="spearman")                          ###
cor.test(mates$mean_rwl.1,mates$belly.avg.bright,method="spearman")                 
cor.test(mates$mean_rwl.1,mates$breast.avg.bright,method="spearman")                
cor.test(mates$mean_rwl.1,mates$throat.avg.bright,method="spearman")
cor.test(mates$mean_rwl.1,mates$vent.avg.bright,method="spearman")                  
cor.test(mates$mean_rwl.1,mates$bird_mass,method="spearman")


# FEMALE variables against MALE right tail streamer

par(mfrow=c(4,3))

plot(mates$mean_rts,mates$proportion.rustica.ancestry.1,pch=20,col=mates$location)
plot(mates$mean_rts,mates$d2h.1,pch=20,col=mates$location)
plot(mates$mean_rts,mates$d13C_VPDB.1,pch=20,col=mates$location)
plot(mates$mean_rts,mates$mean_rwl.1,pch=20,col=mates$location)
plot(mates$mean_rts,mates$mean_lts.1,pch=20,col=mates$location)
plot(mates$mean_rts,mates$belly.avg.bright.1,pch=20,col=mates$location)
plot(mates$mean_rts,mates$breast.avg.bright.1,pch=20,col=mates$location)
plot(mates$mean_rts,mates$throat.avg.bright.1,pch=20,col=mates$location)
plot(mates$mean_rts,mates$vent.avg.bright.1,pch=20,col=mates$location)
plot(mates$mean_rts,mates$bird_mass.1,pch=20,col=mates$location)

cor.test(mates$mean_rts,mates$proportion.rustica.ancestry.1,method="spearman")       ###
cor.test(mates$mean_rts,mates$d2h.1,method="spearman")                         
cor.test(mates$mean_rts,mates$d13C_VPDB.1,method="spearman")                                      
cor.test(mates$mean_rts,mates$mean_rwl.1,method="spearman")                          ###
cor.test(mates$mean_rts,mates$mean_lts.1,method="spearman")                          ###
cor.test(mates$mean_rts,mates$belly.avg.bright.1,method="spearman")                  
cor.test(mates$mean_rts,mates$breast.avg.bright.1,method="spearman")                 
cor.test(mates$mean_rts,mates$throat.avg.bright.1,method="spearman")                 
cor.test(mates$mean_rts,mates$vent.avg.bright.1,method="spearman")                         
cor.test(mates$mean_rts,mates$bird_mass.1,method="spearman")                         ###

# MALE variables against FEMALE right tail streamer

par(mfrow=c(4,3))

plot(mates$mean_rts.1,mates$proportion.rustica.ancestry,pch=20,col=mates$location)
plot(mates$mean_rts.1,mates$d2h,pch=20,col=mates$location)
plot(mates$mean_rts.1,mates$d13C_VPDB,pch=20,col=mates$location)
plot(mates$mean_rts.1,mates$mean_rwl,pch=20,col=mates$location)
plot(mates$mean_rts.1,mates$mean_lts,pch=20,col=mates$location)
plot(mates$mean_rts.1,mates$belly.avg.bright,pch=20,col=mates$location)
plot(mates$mean_rts.1,mates$breast.avg.bright,pch=20,col=mates$location)
plot(mates$mean_rts.1,mates$throat.avg.bright,pch=20,col=mates$location)
plot(mates$mean_rts.1,mates$vent.avg.bright,pch=20,col=mates$location)
plot(mates$mean_rts.1,mates$bird_mass,pch=20,col=mates$location)

cor.test(mates$mean_rts.1,mates$proportion.rustica.ancestry,method="spearman")       
cor.test(mates$mean_rts.1,mates$d2h,method="spearman")                                             
cor.test(mates$mean_rts.1,mates$d13C_VPDB,method="spearman")                                    
cor.test(mates$mean_rts.1,mates$mean_rwl,method="spearman")                          
cor.test(mates$mean_rts.1,mates$mean_lts,method="spearman")                          ###
cor.test(mates$mean_rts.1,mates$belly.avg.bright,method="spearman")                 
cor.test(mates$mean_rts.1,mates$breast.avg.bright,method="spearman")                
cor.test(mates$mean_rts.1,mates$throat.avg.bright,method="spearman")
cor.test(mates$mean_rts.1,mates$vent.avg.bright,method="spearman")                  
cor.test(mates$mean_rts.1,mates$bird_mass,method="spearman")                         ###


# FEMALE variables against MALE belly average brightness

par(mfrow=c(4,3))

plot(mates$belly.avg.bright,mates$proportion.rustica.ancestry.1,pch=20,col=mates$location)
plot(mates$belly.avg.bright,mates$d2h.1,pch=20,col=mates$location)
plot(mates$belly.avg.bright,mates$d13C_VPDB.1,pch=20,col=mates$location)
plot(mates$belly.avg.bright,mates$mean_rwl.1,pch=20,col=mates$location)
plot(mates$belly.avg.bright,mates$mean_lts.1,pch=20,col=mates$location)
plot(mates$belly.avg.bright,mates$mean_rts.1,pch=20,col=mates$location)
plot(mates$belly.avg.bright,mates$breast.avg.bright.1,pch=20,col=mates$location)
plot(mates$belly.avg.bright,mates$throat.avg.bright.1,pch=20,col=mates$location)
plot(mates$belly.avg.bright,mates$vent.avg.bright.1,pch=20,col=mates$location)
plot(mates$belly.avg.bright,mates$bird_mass.1,pch=20,col=mates$location)

cor.test(mates$belly.avg.bright,mates$proportion.rustica.ancestry.1,method="spearman")       
cor.test(mates$belly.avg.bright,mates$d2h.1,method="spearman")                         
cor.test(mates$belly.avg.bright,mates$d13C_VPDB.1,method="spearman")                                      
cor.test(mates$belly.avg.bright,mates$mean_rwl.1,method="spearman")                          
cor.test(mates$belly.avg.bright,mates$mean_lts.1,method="spearman")                          ###
cor.test(mates$belly.avg.bright,mates$mean_rts.1,method="spearman")                  
cor.test(mates$belly.avg.bright,mates$breast.avg.bright.1,method="spearman")                 
cor.test(mates$belly.avg.bright,mates$throat.avg.bright.1,method="spearman")                 ###               
cor.test(mates$belly.avg.bright,mates$vent.avg.bright.1,method="spearman")                         
cor.test(mates$belly.avg.bright,mates$bird_mass.1,method="spearman")                         

# MALE variables against FEMALE belly average brightness

par(mfrow=c(4,3))

plot(mates$belly.avg.bright.1,mates$proportion.rustica.ancestry,pch=20,col=mates$location)
plot(mates$belly.avg.bright.1,mates$d2h,pch=20,col=mates$location)
plot(mates$belly.avg.bright.1,mates$d13C_VPDB,pch=20,col=mates$location)
plot(mates$belly.avg.bright.1,mates$mean_rwl,pch=20,col=mates$location)
plot(mates$belly.avg.bright.1,mates$mean_lts,pch=20,col=mates$location)
plot(mates$belly.avg.bright.1,mates$mean_rts,pch=20,col=mates$location)
plot(mates$belly.avg.bright.1,mates$breast.avg.bright,pch=20,col=mates$location)
plot(mates$belly.avg.bright.1,mates$throat.avg.bright,pch=20,col=mates$location)
plot(mates$belly.avg.bright.1,mates$vent.avg.bright,pch=20,col=mates$location)
plot(mates$belly.avg.bright.1,mates$bird_mass,pch=20,col=mates$location)

cor.test(mates$belly.avg.bright.1,mates$proportion.rustica.ancestry,method="spearman")       
cor.test(mates$belly.avg.bright.1,mates$d2h,method="spearman")                               ###                                               
cor.test(mates$belly.avg.bright.1,mates$d13C_VPDB,method="spearman")                                    
cor.test(mates$belly.avg.bright.1,mates$mean_rwl,method="spearman")                          
cor.test(mates$belly.avg.bright.1,mates$mean_lts,method="spearman")                          
cor.test(mates$belly.avg.bright.1,mates$mean_rts,method="spearman")                 
cor.test(mates$belly.avg.bright.1,mates$breast.avg.bright,method="spearman")                
cor.test(mates$belly.avg.bright.1,mates$throat.avg.bright,method="spearman")
cor.test(mates$belly.avg.bright.1,mates$vent.avg.bright,method="spearman")                  
cor.test(mates$belly.avg.bright.1,mates$bird_mass,method="spearman")                         ###


# FEMALE variables against MALE throat average brightness

par(mfrow=c(4,3))

plot(mates$throat.avg.bright,mates$proportion.rustica.ancestry.1,pch=20,col=mates$location)
plot(mates$throat.avg.bright,mates$d2h.1,pch=20,col=mates$location)
plot(mates$throat.avg.bright,mates$d13C_VPDB.1,pch=20,col=mates$location)
plot(mates$throat.avg.bright,mates$mean_rwl.1,pch=20,col=mates$location)
plot(mates$throat.avg.bright,mates$mean_lts.1,pch=20,col=mates$location)
plot(mates$throat.avg.bright,mates$mean_rts.1,pch=20,col=mates$location)
plot(mates$throat.avg.bright,mates$belly.avg.bright.1,pch=20,col=mates$location)
plot(mates$throat.avg.bright,mates$breast.avg.bright.1,pch=20,col=mates$location)
plot(mates$throat.avg.bright,mates$vent.avg.bright.1,pch=20,col=mates$location)
plot(mates$throat.avg.bright,mates$bird_mass.1,pch=20,col=mates$location)

cor.test(mates$throat.avg.bright,mates$proportion.rustica.ancestry.1,method="spearman")       
cor.test(mates$throat.avg.bright,mates$d2h.1,method="spearman")                         
cor.test(mates$throat.avg.bright,mates$d13C_VPDB.1,method="spearman")                                      
cor.test(mates$throat.avg.bright,mates$mean_rwl.1,method="spearman")                          
cor.test(mates$throat.avg.bright,mates$mean_lts.1,method="spearman")                        ###
cor.test(mates$throat.avg.bright,mates$mean_rts.1,method="spearman")      
cor.test(mates$throat.avg.bright,mates$belly.avg.bright.1,method="spearman")                              
cor.test(mates$throat.avg.bright,mates$breast.avg.bright.1,method="spearman")                 
cor.test(mates$throat.avg.bright,mates$vent.avg.bright.1,method="spearman")                         
cor.test(mates$throat.avg.bright,mates$bird_mass.1,method="spearman")                       ###           

# MALE variables against FEMALE throat average brightness

par(mfrow=c(4,3))

plot(mates$throat.avg.bright.1,mates$proportion.rustica.ancestry,pch=20,col=mates$location)
plot(mates$throat.avg.bright.1,mates$d2h,pch=20,col=mates$location)
plot(mates$throat.avg.bright.1,mates$d13C_VPDB,pch=20,col=mates$location)
plot(mates$throat.avg.bright.1,mates$mean_rwl,pch=20,col=mates$location)
plot(mates$throat.avg.bright.1,mates$mean_lts,pch=20,col=mates$location)
plot(mates$throat.avg.bright.1,mates$mean_rts,pch=20,col=mates$location)
plot(mates$throat.avg.bright.1,mates$belly.avg.bright,pch=20,col=mates$location)
plot(mates$throat.avg.bright.1,mates$breast.avg.bright,pch=20,col=mates$location)
plot(mates$throat.avg.bright.1,mates$vent.avg.bright,pch=20,col=mates$location)
plot(mates$throat.avg.bright.1,mates$bird_mass,pch=20,col=mates$location)

cor.test(mates$throat.avg.bright.1,mates$proportion.rustica.ancestry,method="spearman")       
cor.test(mates$throat.avg.bright.1,mates$d2h,method="spearman")                                                                            
cor.test(mates$throat.avg.bright.1,mates$d13C_VPDB,method="spearman")                         ###                                 
cor.test(mates$throat.avg.bright.1,mates$mean_rwl,method="spearman")                          
cor.test(mates$throat.avg.bright.1,mates$mean_lts,method="spearman")                          
cor.test(mates$throat.avg.bright.1,mates$mean_rts,method="spearman")   
cor.test(mates$throat.avg.bright.1,mates$belly.avg.bright,method="spearman")                  ###
cor.test(mates$throat.avg.bright.1,mates$breast.avg.bright,method="spearman")                
cor.test(mates$throat.avg.bright.1,mates$vent.avg.bright,method="spearman")                   ###         
cor.test(mates$throat.avg.bright.1,mates$bird_mass,method="spearman")                         






### Evaluate relationships for specific localities--------------------------

lan <- mates[which(mates$location=='Lanzhou'),]
zha <- mates[which(mates$location=='Zhangye'),]
jiu <- mates[which(mates$location=='Jiuquan'),]

## Lanzhou
## Isotopes

par(mfrow=c(1,2))
plot(lan$d2h,lan$d2h.1,pch=20)
plot(lan$d13C_VPDB,lan$d13C_VPDB.1,pch=20)

# Are isotope variables correlated for lan?
cor.test(lan$d2h,lan$d2h.1)
# No.
cor.test(lan$d13C_VPDB,lan$d13C_VPDB.1)
# No.

## Ancestry

par(mfrow=c(1,2))
plot(lan$proportion.rustica.ancestry,lan$proportion.rustica.ancestry.1,pch=20)
plot(lan$proportion.gutturalis.ancestry,lan$proportion.gutturalis.ancestry.1,pch=20)

# Are ancestry coefficients correlated?
cor.test(lan$proportion.rustica.ancestry,lan$proportion.rustica.ancestry.1)
# No.
cor.test(lan$proportion.gutturalis.ancestry,lan$proportion.gutturalis.ancestry.1)
# No.

## Wing/tail streamer lengths

par(mfrow=c(1,3))
plot(lan$mean_rwl,lan$mean_rwl.1,pch=20)
plot(lan$mean_rts,lan$mean_rts.1,pch=20)
plot(lan$mean_lts,lan$mean_lts.1,pch=20)

# Are wing/tail streamer lengths correlated?
cor.test(lan$mean_rwl,lan$mean_rwl.1)
# No.
cor.test(lan$mean_rts,lan$mean_rts.1)
# Yes, negative.
cor.test(lan$mean_lts,lan$mean_lts.1)
# Yes, negative.

## Belly feathers

par(mfrow=c(1,4))
plot(lan$belly.total.bright,lan$belly.total.bright.1,pch=20)
plot(lan$belly.avg.bright,lan$belly.avg.bright.1,pch=20)
plot(lan$belly.hue,lan$belly.hue.1,pch=20)
plot(lan$belly.chrom,lan$belly.chrom.1,pch=20)

# Are belly feather variables correlated?
cor.test(lan$belly.total.bright,lan$belly.total.bright.1)
# No.
cor.test(lan$belly.avg.bright,lan$belly.avg.bright.1)
# No.
cor.test(lan$belly.hue,lan$belly.hue.1)
# Yes.
cor.test(lan$belly.chrom,lan$belly.chrom.1)
# No.

## Throat feathers

par(mfrow=c(1,4))
plot(lan$throat.total.bright,lan$throat.total.bright.1,pch=20)
plot(lan$throat.avg.bright,lan$throat.avg.bright.1,pch=20)
plot(lan$throat.hue,lan$throat.hue.1,pch=20)
plot(lan$throat.chrom,lan$throat.chrom.1,pch=20)

# Are throat feather variables correlated?
cor.test(lan$throat.total.bright,lan$throat.total.bright.1)
# No.
cor.test(lan$throat.avg.bright,lan$throat.avg.bright.1)
# No.
cor.test(lan$throat.hue,lan$throat.hue.1)
# No.
cor.test(lan$throat.chrom,lan$throat.chrom.1)
# No.

## Breast feathers

par(mfrow=c(1,4))
plot(lan$breast.total.bright,lan$breast.total.bright.1,pch=20)
plot(lan$breast.avg.bright,lan$breast.avg.bright.1,pch=20)
plot(lan$breast.hue,lan$breast.hue.1,pch=20)
plot(lan$breast.chrom,lan$breast.chrom.1,pch=20)

# Are breast feather variables correlated?
cor.test(lan$breast.total.bright,lan$breast.total.bright.1)
# No.
cor.test(lan$breast.avg.bright,lan$breast.avg.bright.1)
# No.
cor.test(lan$breast.hue,lan$breast.hue.1)
# No.
cor.test(lan$breast.chrom,lan$breast.chrom.1)
# No.

## Vent feathers

par(mfrow=c(1,4))
plot(lan$vent.total.bright,lan$vent.total.bright.1,pch=20)
plot(lan$vent.avg.bright,lan$vent.avg.bright.1,pch=20)
plot(lan$vent.hue,lan$vent.hue.1,pch=20)
plot(lan$vent.chrom,lan$vent.chrom.1,pch=20)

# Are vent feather variables correlated?
cor.test(lan$vent.total.bright,lan$vent.total.bright.1)
# No.
cor.test(lan$vent.avg.bright,lan$vent.avg.bright.1)
# No.
cor.test(lan$vent.hue,lan$vent.hue.1)
# Yes.
cor.test(lan$vent.chrom,lan$vent.chrom.1)
# No.

## Bird mass

par(mfrow=c(1,1))
plot(lan$bird_mass,lan$bird_mass.1,pch=20)

# Is bird mass correlated?
cor.test(lan$bird_mass,lan$bird_mass.1)
# No.


## Zhangye---
## Isotopes

par(mfrow=c(1,2))
plot(zha$d2h,zha$d2h.1,pch=20)
plot(zha$d13C_VPDB,zha$d13C_VPDB.1,pch=20)

# Are isotope variables correlated for zha?
cor.test(zha$d2h,zha$d2h.1)
# No.
cor.test(zha$d13C_VPDB,zha$d13C_VPDB.1)
# No.

## Ancestry

par(mfrow=c(1,2))
plot(zha$proportion.rustica.ancestry,zha$proportion.rustica.ancestry.1,pch=20)
plot(zha$proportion.gutturalis.ancestry,zha$proportion.gutturalis.ancestry.1,pch=20)

# Are ancestry coefficients correlated?
cor.test(zha$proportion.rustica.ancestry,zha$proportion.rustica.ancestry.1)
# No.
cor.test(zha$proportion.gutturalis.ancestry,zha$proportion.gutturalis.ancestry.1)
# No.

## Wing/tail streamer lengths

par(mfrow=c(1,3))
plot(zha$mean_rwl,zha$mean_rwl.1,pch=20)
plot(zha$mean_rts,zha$mean_rts.1,pch=20)
plot(zha$mean_lts,zha$mean_lts.1,pch=20)

# Are wing/tail streamer lengths correlated?
cor.test(zha$mean_rwl,zha$mean_rwl.1)
# No.
cor.test(zha$mean_rts,zha$mean_rts.1)
# Yes, negative.
cor.test(zha$mean_lts,zha$mean_lts.1)
# No.

## Belly feathers

par(mfrow=c(1,4))
plot(zha$belly.total.bright,zha$belly.total.bright.1,pch=20)
plot(zha$belly.avg.bright,zha$belly.avg.bright.1,pch=20)
plot(zha$belly.hue,zha$belly.hue.1,pch=20)
plot(zha$belly.chrom,zha$belly.chrom.1,pch=20)

# Are belly feather variables correlated?
cor.test(zha$belly.total.bright,zha$belly.total.bright.1)
# No.
cor.test(zha$belly.avg.bright,zha$belly.avg.bright.1)
# No.
cor.test(zha$belly.hue,zha$belly.hue.1)
# Yes.
cor.test(zha$belly.chrom,zha$belly.chrom.1)
# No.

## Throat feathers

par(mfrow=c(1,4))
plot(zha$throat.total.bright,zha$throat.total.bright.1,pch=20)
plot(zha$throat.avg.bright,zha$throat.avg.bright.1,pch=20)
plot(zha$throat.hue,zha$throat.hue.1,pch=20)
plot(zha$throat.chrom,zha$throat.chrom.1,pch=20)

# Are throat feather variables correlated?
cor.test(zha$throat.total.bright,zha$throat.total.bright.1)
# No.
cor.test(zha$throat.avg.bright,zha$throat.avg.bright.1)
# No.
cor.test(zha$throat.hue,zha$throat.hue.1)
# No.
cor.test(zha$throat.chrom,zha$throat.chrom.1)
# No.

## Breast feathers

par(mfrow=c(1,4))
plot(zha$breast.total.bright,zha$breast.total.bright.1,pch=20)
plot(zha$breast.avg.bright,zha$breast.avg.bright.1,pch=20)
plot(zha$breast.hue,zha$breast.hue.1,pch=20)
plot(zha$breast.chrom,zha$breast.chrom.1,pch=20)

# Are breast feather variables correlated?
cor.test(zha$breast.total.bright,zha$breast.total.bright.1)
# No.
cor.test(zha$breast.avg.bright,zha$breast.avg.bright.1)
# No.
cor.test(zha$breast.hue,zha$breast.hue.1)
# No.
cor.test(zha$breast.chrom,zha$breast.chrom.1)
# No.

## Vent feathers

par(mfrow=c(1,4))
plot(zha$vent.total.bright,zha$vent.total.bright.1,pch=20)
plot(zha$vent.avg.bright,zha$vent.avg.bright.1,pch=20)
plot(zha$vent.hue,zha$vent.hue.1,pch=20)
plot(zha$vent.chrom,zha$vent.chrom.1,pch=20)

# Are vent feather variables correlated?
cor.test(zha$vent.total.bright,zha$vent.total.bright.1)
# No.
cor.test(zha$vent.avg.bright,zha$vent.avg.bright.1)
# No.
cor.test(zha$vent.hue,zha$vent.hue.1)
# No.
cor.test(zha$vent.chrom,zha$vent.chrom.1)
# No.

## Bird mass

par(mfrow=c(1,1))
plot(zha$bird_mass,zha$bird_mass.1,pch=20)

# Is bird mass correlated?
cor.test(zha$bird_mass,zha$bird_mass.1)
# No.


### Evaluate relationships by year----------------------------

y16 <- mates[which(mates$year=='2016'),]
y17 <- mates[which(mates$year=='2017'),]

## 2016---
## Isotopes

par(mfrow=c(1,2))
plot(y16$d2h,y16$d2h.1,pch=20)
plot(y16$d13C_VPDB,y16$d13C_VPDB.1,pch=20)

# Are isotope variables correlated for y16?
cor.test(y16$d2h,y16$d2h.1)
# No.
cor.test(y16$d13C_VPDB,y16$d13C_VPDB.1)
# No.

## Ancestry

par(mfrow=c(1,2))
plot(y16$proportion.rustica.ancestry,y16$proportion.rustica.ancestry.1,pch=20)
plot(y16$proportion.gutturalis.ancestry,y16$proportion.gutturalis.ancestry.1,pch=20)

# Are ancestry coefficients correlated?
cor.test(y16$proportion.rustica.ancestry,y16$proportion.rustica.ancestry.1)
# No.
cor.test(y16$proportion.gutturalis.ancestry,y16$proportion.gutturalis.ancestry.1)
# No.

## Wing/tail streamer lengths

par(mfrow=c(1,3))
plot(y16$mean_rwl,y16$mean_rwl.1,pch=20)
plot(y16$mean_rts,y16$mean_rts.1,pch=20)
plot(y16$mean_lts,y16$mean_lts.1,pch=20)

# Are wing/tail streamer lengths correlated?
cor.test(y16$mean_rwl,y16$mean_rwl.1)
# No.
cor.test(y16$mean_rts,y16$mean_rts.1)
# Yes, negative.
cor.test(y16$mean_lts,y16$mean_lts.1)
# No.

## Belly feathers

par(mfrow=c(1,4))
plot(y16$belly.total.bright,y16$belly.total.bright.1,pch=20)
plot(y16$belly.avg.bright,y16$belly.avg.bright.1,pch=20)
plot(y16$belly.hue,y16$belly.hue.1,pch=20)
plot(y16$belly.chrom,y16$belly.chrom.1,pch=20)

# Are belly feather variables correlated?
cor.test(y16$belly.total.bright,y16$belly.total.bright.1)
# No.
cor.test(y16$belly.avg.bright,y16$belly.avg.bright.1)
# No.
cor.test(y16$belly.hue,y16$belly.hue.1)
# No.
cor.test(y16$belly.chrom,y16$belly.chrom.1)
# Yes.

## Throat feathers

par(mfrow=c(1,4))
plot(y16$throat.total.bright,y16$throat.total.bright.1,pch=20)
plot(y16$throat.avg.bright,y16$throat.avg.bright.1,pch=20)
plot(y16$throat.hue,y16$throat.hue.1,pch=20)
plot(y16$throat.chrom,y16$throat.chrom.1,pch=20)

# Are throat feather variables correlated?
cor.test(y16$throat.total.bright,y16$throat.total.bright.1)
# No.
cor.test(y16$throat.avg.bright,y16$throat.avg.bright.1)
# No.
cor.test(y16$throat.hue,y16$throat.hue.1)
# No.
cor.test(y16$throat.chrom,y16$throat.chrom.1)
# No.

## Breast feathers

par(mfrow=c(1,4))
plot(y16$breast.total.bright,y16$breast.total.bright.1,pch=20)
plot(y16$breast.avg.bright,y16$breast.avg.bright.1,pch=20)
plot(y16$breast.hue,y16$breast.hue.1,pch=20)
plot(y16$breast.chrom,y16$breast.chrom.1,pch=20)

# Are breast feather variables correlated?
cor.test(y16$breast.total.bright,y16$breast.total.bright.1)
# No.
cor.test(y16$breast.avg.bright,y16$breast.avg.bright.1)
# No.
cor.test(y16$breast.hue,y16$breast.hue.1)
# No.
cor.test(y16$breast.chrom,y16$breast.chrom.1)
# Yes.

## Vent feathers

par(mfrow=c(1,4))
plot(y16$vent.total.bright,y16$vent.total.bright.1,pch=20)
plot(y16$vent.avg.bright,y16$vent.avg.bright.1,pch=20)
plot(y16$vent.hue,y16$vent.hue.1,pch=20)
plot(y16$vent.chrom,y16$vent.chrom.1,pch=20)

# Are vent feather variables correlated?
cor.test(y16$vent.total.bright,y16$vent.total.bright.1)
# No.
cor.test(y16$vent.avg.bright,y16$vent.avg.bright.1)
# No.
cor.test(y16$vent.hue,y16$vent.hue.1)
# No.
cor.test(y16$vent.chrom,y16$vent.chrom.1)
# No.

## Bird mass

par(mfrow=c(1,1))
plot(y16$bird_mass,y16$bird_mass.1,pch=20)

# Is bird mass correlated?
cor.test(y16$bird_mass,y16$bird_mass.1)
# Yes.


## 2017---
## Isotopes

par(mfrow=c(1,2))
plot(y17$d2h,y17$d2h.1,pch=20)
plot(y17$d13C_VPDB,y17$d13C_VPDB.1,pch=20)

# Are isotope variables correlated for y17?
cor.test(y17$d2h,y17$d2h.1)
# No.
cor.test(y17$d13C_VPDB,y17$d13C_VPDB.1)
# No.

## Ancestry

par(mfrow=c(1,2))
plot(y17$proportion.rustica.ancestry,y17$proportion.rustica.ancestry.1,pch=20)
plot(y17$proportion.gutturalis.ancestry,y17$proportion.gutturalis.ancestry.1,pch=20)

# Are ancestry coefficients correlated?
cor.test(y17$proportion.rustica.ancestry,y17$proportion.rustica.ancestry.1)
# No.
cor.test(y17$proportion.gutturalis.ancestry,y17$proportion.gutturalis.ancestry.1)
# No.

## Wing/tail streamer lengths

par(mfrow=c(1,3))
plot(y17$mean_rwl,y17$mean_rwl.1,pch=20)
plot(y17$mean_rts,y17$mean_rts.1,pch=20)
plot(y17$mean_lts,y17$mean_lts.1,pch=20)

# Are wing/tail streamer lengths correlated?
cor.test(y17$mean_rwl,y17$mean_rwl.1)
# No.
cor.test(y17$mean_rts,y17$mean_rts.1)
# Yes, negative.
cor.test(y17$mean_lts,y17$mean_lts.1)
# No.

## Belly feathers

par(mfrow=c(1,4))
plot(y17$belly.total.bright,y17$belly.total.bright.1,pch=20)
plot(y17$belly.avg.bright,y17$belly.avg.bright.1,pch=20)
plot(y17$belly.hue,y17$belly.hue.1,pch=20)
plot(y17$belly.chrom,y17$belly.chrom.1,pch=20)

# Are belly feather variables correlated?
cor.test(y17$belly.total.bright,y17$belly.total.bright.1)
# No.
cor.test(y17$belly.avg.bright,y17$belly.avg.bright.1)
# No.
cor.test(y17$belly.hue,y17$belly.hue.1)
# No.
cor.test(y17$belly.chrom,y17$belly.chrom.1)
# Yes.

## Throat feathers

par(mfrow=c(1,4))
plot(y17$throat.total.bright,y17$throat.total.bright.1,pch=20)
plot(y17$throat.avg.bright,y17$throat.avg.bright.1,pch=20)
plot(y17$throat.hue,y17$throat.hue.1,pch=20)
plot(y17$throat.chrom,y17$throat.chrom.1,pch=20)

# Are throat feather variables correlated?
cor.test(y17$throat.total.bright,y17$throat.total.bright.1)
# No.
cor.test(y17$throat.avg.bright,y17$throat.avg.bright.1)
# No.
cor.test(y17$throat.hue,y17$throat.hue.1)
# No.
cor.test(y17$throat.chrom,y17$throat.chrom.1)
# No.

## Breast feathers

par(mfrow=c(1,4))
plot(y17$breast.total.bright,y17$breast.total.bright.1,pch=20)
plot(y17$breast.avg.bright,y17$breast.avg.bright.1,pch=20)
plot(y17$breast.hue,y17$breast.hue.1,pch=20)
plot(y17$breast.chrom,y17$breast.chrom.1,pch=20)

# Are breast feather variables correlated?
cor.test(y17$breast.total.bright,y17$breast.total.bright.1)
# No.
cor.test(y17$breast.avg.bright,y17$breast.avg.bright.1)
# No.
cor.test(y17$breast.hue,y17$breast.hue.1)
# No.
cor.test(y17$breast.chrom,y17$breast.chrom.1)
# No.

## Vent feathers

par(mfrow=c(1,4))
plot(y17$vent.total.bright,y17$vent.total.bright.1,pch=20)
plot(y17$vent.avg.bright,y17$vent.avg.bright.1,pch=20)
plot(y17$vent.hue,y17$vent.hue.1,pch=20)
plot(y17$vent.chrom,y17$vent.chrom.1,pch=20)

# Are vent feather variables correlated?
cor.test(y17$vent.total.bright,y17$vent.total.bright.1)
# No.
cor.test(y17$vent.avg.bright,y17$vent.avg.bright.1)
# No.
cor.test(y17$vent.hue,y17$vent.hue.1)
# No.
cor.test(y17$vent.chrom,y17$vent.chrom.1)
# No.

## Bird mass

par(mfrow=c(1,1))
plot(y17$bird_mass,y17$bird_mass.1,pch=20)

# Is bird mass correlated?
cor.test(y17$bird_mass,y17$bird_mass.1)
# No.
