# calculating nucleotide diversity

library(radiator)

data<-tidy_vcf("olfusa.filtered.vcf", strata="strata.olfusa.clean.tsv")
pi(data)


PI<-read.table("pi.individuals.tsv", header = T)

PI$POP_ID = factor(PI$POP_ID, c("HIN", "MID", "FRE", "VAR", "HES", "SOG", "OLF",
                                "HVI", "EFR", "FUS", "OXA", "THI", "ULF", "VIL",
                                "THV", "LEI", "SLP"))

my.colors<-c("#a50026", "#a50026", "#a50026", "#a50026", "#f46d43",
                      "#f46d43", "#f46d43", "#f46d43", "#fee090", "#fee090", 
                      "#fee090", "#fee090", "#fee090", "#fee090", "#fee090", 

                                            "#74add1", "#313695")
par(mar=c(5,5,2,2))
par(family="Times")
boxplot(PI$PI~PI$POP_ID,
        ylim=c(0,13),
        xlab="Location",
        ylab=expression(paste("Nucleotide diversity (π) x 10"^"-4")),
        col=my.colors,
        density=c(5,10,20,30,7,
                  3,3,3,3,3,
                  3,3,3,3,3,3,3)
)