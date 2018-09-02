################################
## READ AND PARSE RECON MODEL ##
################################

## recon is assumed to be available by the path stored in the model.file variable

recon           <- read.table(model.file, sep="\t")
colnames(recon) <- c("name", "reaction", "min", "max", "rules", "note")
rownames(recon) <- paste("i", 1:nrow(recon), sep="")

cat(paste("Number of all reactions:\t", nrow(recon),"\n"))

## ENZYMATIC REACTIONS SELECTION (HAVING A GENETIC RULE)
enzymatic.reactions <- recon[recon$rules != "",]
cat(paste("Number of enzymatic reactions:\t", nrow(enzymatic.reactions),"\n"))

## STATISTICS ON GENE HGNC IDs
genes.by.rules      <- strsplit(as.character(recon$rules), "[^0-9]+")
genes.per.reaction  <- sapply(genes.by.rules, length)

all.genes.in.rules  <- na.omit(as.numeric(unlist(genes.by.rules)))
occurences.of.gene  <- table(all.genes.in.rules)
#sort(occurences.of.gene, decreasing = T)[1:15]

unique.genes        <- unique(all.genes.in.rules)
cat(paste("Number of unique genes\namong enzymatic reactions:\t", length(unique.genes), "\n"))
