###################
# Author: Qunhua Li
# Email: qli@stat.berkeley.edu
# Affiliation: Peter Bickel's group, Dept of Statistics, UC Berkeley
###################
# 12-06-2010
# Move mixing proportion p to command line
#
# 10-10-2010
#
# This program runs IDR analysis for matched peaks
# It can handle ties and discreteness in the data
#  
# NOTE: this is still a beta version, not tested extensively yet
#       Comments and feedbacks are welcome
#
# Command:
#
# Rscript batch-IDR-mappedInput-discrete.r matched.file outfile.prefix include.all bin.size p
#
# Another option is to hardcode the fields above in this file and use the following command (It can solve the problem of failing to find package adapt)
#
# R CMD BATCH batch-IDR-mappedInput-discrete.r 
#
# Input:
#  matched.file: filename of matched peaks
#      The format of input file is "ID1 score1 ID2  score2"
#      where each row represents a matched pair of peaks, with
#      ID1=-1 and score1=-1 if a peak is only on replicate2 but not on
#      replicate1 and vice versa. 
#  outfile.prefix: prefix of output file
#  include.all: T: include all peaks; F: only include overlapping ones 
#  bin.size: controls the bins to take integrations on,
#            If bin.size=n, integrate all bins with num of observations>n and
#            treat obs in bins with num of obs<=n as individual obs
#            The larger bin.size is, the coarser the results will be, but
#            the fast it is. The main impact is on 
#            When n=1, get exact results, but super slow for large datasets
#            (can take days to run)
#            Suggested value: 5 - 20
#  p: starting values for proportion of peaks that are reproducible
#
#  Output: Let *=outfile.prefix
#  *-Rout.txt : excution log
#  *-uri.sav: output R object for plotting Psi and Psi' curve,
#             need use load() to R to read
#  *-em.sav: Output R object from copula mixture model,
#            need use load() to load to R
#  *-overlapped-peaks.txt: peaks appear on both replicates,
#            with local idr score for each pair of matched peaks 
#  *-npeaks-aboveIDR.txt: Number of peaks selected at IDR cutoff in the
#     range of (0.01-0.25), also the corresponding significance threshold
#     on each replicate at each IDR cutoff
#  *-plot.ps: IDR curve and a diagnosis plot to show the distribution of
#             local idr
#             In the diagnosis plot, grey means local idr is close to 0
#               black means local idr is close to 1


args <- commandArgs(trailingOnly=T)

matched.file <- args[1]
output.prefix <- args[2]
include.all <- as.logical(args[3]) # T: include all peaks; F: only include overlapping ones
bin.size <- as.numeric(args[4]) # controls the bins to take integrations on,
p <- as.numeric(args[5]) # proportion of reproducible peaks

source("functions-09-12-2010.r",local=TRUE)

#matched.file <- "~/ENCODE/RNA-discrete/data/alex/SJforIDR.txt"
#output.prefix <- "~/ENCODE/RNA-discrete/data/alex/result-SJforIDR"
#include.all <- F # include peaks tht are present on only one replicate
#bin.size <- 20

sink(paste(output.prefix, "-Rout.txt", sep=""))

data.file <- read.table(paste(matched.file, sep=""))

x <- data.file[,2]
y <- data.file[,4]

cat(paste(nrow(data.file), "peaks are read.\n"))
cat(paste(sum(x!=-1 & y!=-1), "peaks are present on both replicates.\n\n"))
cat(paste("Include all peaks: ", include.all, "\n"))
cat(paste("Integration is taken for bins with more than", bin.size, "observations.\n"))

if(!include.all){
  is.on.both <- data.file[,2] !=-1 & data.file[,4] != -1
  x <- x[is.on.both]
  y <- y[is.on.both]
}
  
cat(paste("Computing upper rank correlation (URI) for correspondence curves...\n"))

uri.output <- get.uri.matched(x, y)


# save uri output
save(uri.output, file=paste(output.prefix, "-uri.sav", sep=""))

cat(paste("URI is saved at: ", output.prefix, "-uri.sav \n", sep=""))

# EM procedure for inference
cat(paste("Fitting copula mixture model ...\n"))

mu <- 2.6
sigma <- 1.2
rho <- 0.8
#p <- 0.6
eps <- 0.01


em.output <- em.transform.discrete(x, y, mu, sigma, rho, p, eps, bin.size)

save(em.output, file=paste(output.prefix, "-em.sav", sep=""))
#load(paste(output.prefix, "-em.sav", sep=""))
#em.output <- em.output.nomiss

cat(paste("EM is done\n\n"))

#save(em.output, file=paste(output.prefix, "-em-bin5.sav", sep=""))
cat(paste("EM is saved at: ", output.prefix, "-em-bin10.sav \n", sep=""))

# write em output into a file
cat(paste("Estimated parameters in copula mixture model for the following files\n", matched.file))
print(em.output$para)

#for testing
#load(paste(output.prefix, "-em-bin5.sav", sep=""))
#load(paste(output.prefix, "-uri.sav", sep=""))

data.output <- find.idr.IDR(x, y, em.output)
cat(paste("Write peaks and their local idr and IDR to: ", output.prefix, "-pairwise-idr.txt", sep=""))
write.table(data.output$idr.data, file=paste(output.prefix, "-pairwise-idr.txt", sep=""), row.names=F)

cat(paste("Write counts, local idr and IDR of peak categories to: ", output.prefix, "-peaks.txt", sep=""))
write.table(data.output$idr.cat.data, file=paste(output.prefix, "-categories.txt", sep=""), row.names=F)


# number of peaks passing IDR range (0.002-0.2)
cat(paste("Write number of peaks for IDR cutoff 0.002-0.2 to", output.prefix, "-npeaks.txt", sep=""))
IDR.cutoff <- seq(0.002, 0.2, by=0.002)
output.cutoff <- get.npeaks.above.cutoff(data.output$idr.cat.data, IDR.cutoff)
write.table(output.cutoff, file=paste(output.prefix, "-npeaks.txt", sep=""), row.names=F)

sink()

###
### plotting
### 
postscript(paste(output.prefix, "-plot.ps", sep=""))
par(mfrow=c(2,2))
# plot IDR curve
plot.IDR.curve(data.output$idr.cat.data, plotfile=NULL, title.txt="")
#getImage(data.output$idr.cat.data$score1, data.output$idr.cat.data$score2, data.output$idr.cat.data$idr, "idr distribution")

# plot decision boundary
is.dup <- duplicated(cbind(x,y))
index.unique <- which(!is.dup)
is.in <- data.output$idr.cat.data$idr[index.unique] >0.1
rank1 <- rank(data.output$idr.cat.data$score1)[index.unique]
rank2 <- rank(data.output$idr.cat.data$score2)[index.unique]
plot(rank1[is.in], rank2[is.in])
points(rank1[!is.in], rank2[!is.in], col=2)


# plot URI curve
plot.uri.group(list(uri.output$uri.n), NULL, file.name=NULL, 1, title.txt="")
dev.off()

# need put this part into a function and write output to a file





### map idr back to original significance values on the replicates
# Caution if you want to use the mapped original significance values:
# This mapping step uses heuristics that deteriorate the idea of
# reproducibility, because idr and original thresholds are not monotonely
# mapped. It introduces unstability and may keep noise peaks as signal.
# Not recommended, unless you desperately want to have more peaks at the cost
# of irreproducibility!  Interpret the significance values with caution!!!
###
# # You should always use idr as the criterion to choose peaks.
###

######old code############################################

# try different parameters

#mu <- 1.6
#sigma <- 0.5
#rho <- 0.8
#p <- 0.6
#eps <- 0.01

#em.output2 <- em.transform.discrete(x, y, mu, sigma, rho, p, eps, 5)

#save(em.output2, file=paste(output.prefix, "-em2-bin5.sav", sep=""))

# try different parameters
#mu <- 2
#sigma <- 1
#rho <- 0.5
#p <- 0.5
#eps <- 0.01

#em.output3 <- em.transform.discrete(x, y, mu, sigma, rho, p, eps, 5)

#save(em.output3, file=paste(output.prefix, "-em3-bin5.sav", sep=""))

