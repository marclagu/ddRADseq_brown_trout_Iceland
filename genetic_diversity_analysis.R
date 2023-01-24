## Fst calculation

library(assigner)
library(radiator)

data<-tidy_vcf("olfusa.filtered.vcf", strata="strata.olfusa.clean.tsv")

fst.all <- fst_WC84(data, snprelate = FALSE, strata = "strata.olfusa.clean.tsv", 
                    pop.levels = c("EFR", "FUS", "OXA", "THI", "ULF", "VIL", 
                                  "THV", "HIN", "MID", "FRE", "HES", "HVI",
                                  "SOG", "OLF", "VAR", "LEI", "SLP"), pairwise = TRUE, ci = TRUE, iteration.ci = 100)

df <- fst.all$pairwise.fst

write.csv(df, "fst.all")

## He and Ho

library(vcfR)
library(adegenet)
library(hierfstat)
library(diveRsity)

vcf <- read.vcfR("olfusa.filtered.vcf")

x <- vcfR2genind(vcf)

new.pop <-read.table("strata.olfusa.tsv", header=F)
pop(x) <- new.pop$V2

my.colors<-c("#fee090","#a50026", "#fee090", "#f46d43", "#a50026",
                      "#f46d43", "#74add1", "#a50026", "#f46d43", "#fee090", 
                      "#313695", "#f46d43", "#fee090", "#fee090", "#fee090", 
                      "#a50026", "#fee090")

par(mar=c(5.1, 4.1, 4.1, 2.1))

barplot(table(pop(x)), las=3,
        xlab="Individual", ylab="Sample size", col=my.colors)

temp <- summary(x)
str(temp)
barplot(temp$n.by.pop, las=3,
        xlab="Population", ylab="Number of individuals", col=my.colors)

barplot(temp$pop.n.all, las=3,
        xlab="Population", ylab="Number of alleles", col=my.colors)


plot(temp$Hexp, temp$Hobs, pch=20, cex=1, xlim=c(0,0.5), ylim=c(0,0.5), xlab="He", ylab="Ho")
abline(0,1,lty=2)

bartlett.test(list(temp$Hexp, temp$Hobs))

# BARTLETT TEST PER POPULATION
table <- read.csv("Ho_vs_He_per_population.csv", sep = ",", header = TRUE)

bartlett.test(list(table$Ho_EFR, table$He_EFR))
bartlett.test(list(table$Ho_FRE, table$He_FRE))
bartlett.test(list(table$Ho_MID, table$He_MID))
bartlett.test(list(table$Ho_HIN, table$He_HIN))
bartlett.test(list(table$Ho_HES, table$He_HES))
bartlett.test(list(table$Ho_HVI, table$He_HVI))
bartlett.test(list(table$Ho_LEI, table$He_LEI))
bartlett.test(list(table$Ho_OLF, table$He_OLF))
bartlett.test(list(table$Ho_FUS, table$He_FUS))
bartlett.test(list(table$Ho_OXA, table$He_OXA))
bartlett.test(list(table$Ho_SLP, table$He_SLP))
bartlett.test(list(table$Ho_SOG, table$He_SOG))
bartlett.test(list(table$Ho_THI, table$He_THI))
bartlett.test(list(table$Ho_THV, table$He_THV))
bartlett.test(list(table$Ho_ULF, table$He_ULF))
bartlett.test(list(table$Ho_VAR, table$He_VAR))
bartlett.test(list(table$Ho_VIL, table$He_VIL))

# Basic analysis with hierfstat
basic <- basic.stats(x)