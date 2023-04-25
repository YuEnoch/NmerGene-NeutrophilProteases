#' Calculate a quantitative map of teh sequence space of a multiple sequence alignment
#'
#' @param MSA         MSA as a fasta file, matrix of characters, or data frame of characters
#' @param res.prop    Residue property table as a csv, matrix, or data frame (default = default values for charge, hydrophobicity, R-group molecular weight, disorder propensity, occupancy)
#' @param clusters    The number of sequence clusters to test for (default = 1:8)
#' @param clusterPCs  The number of principal components to be used when identifing clusters of sequences (default 5)
#' @param cys         Whether the sequence is cysteine-rich (default = TRUE)
#' @param model      Which shapes of clusters to test for (default = "EII", "VII", "EEI", "VEI", "EVI", "VVI", "EEE", "EVE", "VEE", "VVE", "EEV", "VEV", "EVV", "VVV")
#'
#' @return Generates an object of class "Sequence Alignment Principal Component Analysis" (SAPCA), providing the scales, numericised matrix generated from the MSA. Data is scaled within each property type. The details of the PAC output components are as follows:
#' \itemize{
#'  \item {MSA}             {Input MSA as a matrix of characters}
#'  \item {res.prop}        {Input residue property table as a data frame}
#'  \item {MSA.num.stack}   {Numericised MSA with each property as a separate list item}
#' }
#'
#' @export
#' @examples
#' data(example_MSA)
#' data(residue_properties)
#' SAPCA <- PCA_MSA(example_MSA,residue_properties)
#'
#' @note The `PCA_MSA' function performs multidimensional scaling of a numericised MSA into a reduced space by PCA, to simplify it and aid interpretation.
#'

library("mclust")

#library("rgl")
#library("rglwidget")
#library("ape")
#library("ggplot2")
#library("tidyr")
#library("seqinr")

# Main PCA -------

PCA_MSA <- function(MSA,
                    res.prop,
                    cys        = TRUE,
                    clusterPCs = 1:40,
                    clusters   = 1:15,
                    model      = "VVV",
                    bootstrap  = NULL,
                    jackknife  = NULL,
                    sam        = 0.9,
                    repPCs     = 3,
                    correlation= "pearson"){

  # for some reason, seem to need to explicitly load this library
  library(quietly = TRUE,"mclust")


  ########################.
  # Perform subfunctions #
  ########################.
  MSA       <- read.MSA(MSA = MSA)
  res.prop  <- read.res.prop(res.prop = res.prop,
                             cys      = cys)
  message("generating quantitative sequence space matrix")
  numerical.alignment <- numericise_MSA(MSA      = MSA,
                                        res.prop = res.prop,
                                        cys      = cys)
  message("rotating and projecting sequence space matrix")
  seq.space.PCA <- rotate_seqspace(numerical.alignment = numerical.alignment,
                                   cys                 = cys)
  message("identifying clusters")
    # if no replicates
  if(is.null(bootstrap) && is.null(jackknife)){
    seq.space.clusters <- find_clusters(seq.space.PCA = seq.space.PCA,
                                        clusterPCs    = clusterPCs,
                                        clusters      = clusters,
                                        model         = model)
  }

  #if bootstrap replicates
  if(!is.null(bootstrap)){
    bootstrap.reps <- bootstrap(MSA,
                                res.prop   = res.prop,
                                cys        = cys,
                                clusterPCs = clusterPCs,
                                clusters   = clusters,
                                model      = model,
                                sam        = sam,
                                reps       = bootstrap,
                                PCs        = repPCs,
                                correlation= correlation)
    seq.space.clusters <- find_clusters(seq.space.PCA = seq.space.PCA,
                                        clusterPCs    = clusterPCs,
                                        clusters      = bootstrap.reps$clust.opt,
                                        model         = model)
  }

  #if jackknife replicates
  if(!is.null(jackknife)){
    jackknife.reps <- jackknife(MSA,
                                res.prop   = res.prop,
                                cys        = cys,
                                clusterPCs = clusterPCs,
                                clusters   = clusters,
                                model      = model,
                                sam        = sam,
                                reps       = jackknife,
                                PCs        = repPCs,
                                correlation= correlation)
    seq.space.clusters <- find_clusters(seq.space.PCA = seq.space.PCA,
                                        clusterPCs    = clusterPCs,
                                        clusters      = jackknife.reps$clust.opt,
                                        model         = model)
  }

  #################.
  # Format output #
  #################.
  output <- list(numerical.alignment = numerical.alignment,
                 seq.space.PCA       = seq.space.PCA,
                 seq.space.clusters  = seq.space.clusters,
                 call                = list(MSA        = MSA,
                                            res.prop   = res.prop,
                                            cys        = cys,
                                            clusterPCs = clusterPCs,
                                            clusters   = clusters))

  #if replicates have been done, add them to the output
  if(!is.null(bootstrap)){
    output$bootstrap.reps   <- bootstrap.reps
    output$call$bootstrap   <- bootstrap
    output$call$sam         <- sam
    output$call$repPCs      <- repPCs
    output$call$correlation <- correlation
  }
  if(!is.null(jackknife)){
    output$jackknife.reps   <- jackknife.reps
    output$call$jackknife   <- jackknife
    output$call$sam         <- sam
    output$call$repPCs      <- repPCs
    output$call$correlation <- correlation
  }

  output
}



##################.
# Numericise MSA #
##################.

numericise_MSA <- function(MSA,
                           res.prop,
                           cys){
  seq.names <- rownames(MSA)
  aln.len   <- ncol(MSA)
  res.props <- colnames(res.prop)
  res.avail <- row.names(res.prop)
  # Numericise MSA based on res.prop
  MSA.num.tall <- res.prop[t(MSA),]
  # Name data types
  rownames(MSA.num.tall) <- NULL
  sequence     <- rep(x = seq.names, each  = aln.len)
  residue      <- rep(x = 1:aln.len, times = length(seq.names))
  MSA.num.tall <- cbind(sequence, residue, MSA.num.tall)
  # Stack data into list of matrices
  MSA.num.stack <- NULL
  for (x in 1:length(res.props)) {
    col.names <- paste(1:aln.len,
                       rep(res.props[x],aln.len),
                       sep = ".")
    MSA.num.stack[[res.props[x]]] <- matrix(MSA.num.tall[,x+2],
                                            ncol     = aln.len,
                                            byrow    = TRUE,
                                            dimnames = list(seq.names,
                                                            col.names))
  }
  # Also reflow into single wide matrix
  MSA.num.wide <- MSA.num.stack[[1]]
  for (x in 2:length(res.props)) {
    MSA.num.wide <- cbind(MSA.num.wide, MSA.num.stack[[res.props[x]]])
  }


  ############################.
  # Scaling by property type #
  ############################.

  # Take means and variances of each property type
  prop.means <- NULL
  prop.vars  <- NULL
  for (x in 1:length(res.props)) {
    prop.means[x] <- mean(MSA.num.stack[[x]],na.rm=1)
    prop.vars[x]  <- var(tidyr::gather(data.frame(MSA.num.stack[[x]]))[2],na.rm=1)
  }
  names(prop.means) <- res.props
  names(prop.vars)  <- res.props

  # Scale numericised MSA to prop.means and prop.vars
  MSA.scale.stack <- NULL
  for (x in 1:length(res.props)) {
    MSA.scale.stack[[res.props[x]]] <- (MSA.num.stack[[res.props[x]]]- prop.means[x]) /
      sqrt(prop.vars[x])
  }

  # Replace gaps (currently "NA") with column average
  # Create na.colmean function
  na.colmean<-function(x){
    x[is.na(x)] <- mean(as.matrix(x),na.rm = 1)
    x
  }
  # For each property of MSA.num.stack, apply na.colmean function to each matrix comlumn
  for (x in 1:length(res.props)) {
    MSA.scale.stack[[x]] <- apply(MSA.scale.stack[[x]],2,na.colmean)
  }

  # Also reflow into singe wide matrix for PCA
  MSA.scale.wide <- MSA.scale.stack[[1]]
  for (x in 2:length(res.props)) {
    MSA.scale.wide <- cbind(MSA.scale.wide, MSA.scale.stack[[x]])
  }


  ##################.
  # Alignment list #
  ##################.
  numerical.alignment <- list(MSA             = MSA,
                              res.prop        = res.prop,
                              MSA.num.stack   = MSA.num.stack,
                              MSA.num.wide    = MSA.num.wide,
                              MSA.scale.stack = MSA.scale.stack,
                              MSA.scale.wide  = MSA.scale.wide,
                              prop.means      = prop.means,
                              prop.vars       = prop.vars,
                              seq.names       = seq.names,
                              aln.len         = aln.len)
  numerical.alignment
}


##############################################################.
# Principal Component analysis of numericised sequence space #
##############################################################.

rotate_seqspace <- function(numerical.alignment,
                            cys = TRUE){
  # Prepare data for PCA
  toPCA <- numerical.alignment$MSA.scale.wide
  # If any columns contained gaps only, replace all values with 0
  toPCA[,is.na(colMeans(toPCA))]<-0

  if (cys==FALSE){
    toPCA <- toPCA[,grep("CYS", colnames(toPCA), invert = TRUE)]
  }

  # PCA of data with no extra scaling
  PCA.raw <- stats::prcomp(toPCA)

  seq.space.PCA <- list(coordinates = PCA.raw$x,
                        loadings    = PCA.raw$rotation,
                        centre      = PCA.raw$center,
                        scale       = PCA.raw$scale,
                        stdev       = PCA.raw$sdev)
  seq.space.PCA
}


################################.
# Finding model-based clusters #
################################.

find_clusters <- function(seq.space.PCA,
                          clusterPCs,
                          clusters,
                          model){
  # Mclust calculating multiple cluster numbsers and chooses most likely
  # Different model assume different (non-circular) cluster shapes
  # Different numbers of clusters
  # Manual - http://www.stat.washington.edu/research/reports/2012/tr597.pdf
  # model - http://finzi.psych.upenn.edu/R/library/mclust/html/mclustModelNames.html

  if (length(clusterPCs)==1){
    clusterPCs <- 1:clusterPCs
  }

  # Using Mclust on PCA data
  clusters.raw <- mclust::Mclust(seq.space.PCA$coordinates[,clusterPCs], # which PCs to use to find clusters
                                 prior      = priorControl(),            # starting values for BIC iterations
                                 G          = clusters,                  # number of possible clusters to assess
                                 modelNames = model)                     # model to test

  clusters.min                <- clusters.raw
  clusters.min$classification <- NULL
  clusters.min$G              <- NULL
  clusters.min$z              <- NULL
  clusters.min$call           <- NULL

  seq.space.clusters <- list(classification   = clusters.raw$classification,
                             optimal          = clusters.raw$G,
                             likelihoods      = clusters.raw$z,
                             other            = clusters.min)
  seq.space.clusters
}


################################.
# Bootstral and jackknife reps #
################################.

bootstrap <- function(MSA,
                      res.prop   = res.prop2,
                      cys        = cys,
                      clusterPCs = 1:20,
                      clusters   = 1:15,
                      model      = "VVV",
                      sam        = 0.9,
                      reps       = 10,
                      PCs        = 3,
                      correlation= "pearson"){
  clust.fit      <- NULL
  replicates     <- NULL
  coords         <- NULL
  coords.cor     <- NULL
  coords.av.score<- NULL

  message(paste(reps, "bootstrapped replicates using random",sam,"column samples"))
  for(rep in 1:reps){
    replicates[[rep]] <- PCA_MSA(MSA=MSA[,sample(ncol(MSA), ncol(MSA)*sam)],
                                 res.prop   = res.prop,
                                 clusterPCs = clusterPCs,
                                 clusters   = clusters,
                                 model      = model,
                                 bootstrap  = NULL,
                                 jackknife  = NULL)

    clust.fit<-rbind(clust.fit,as.vector(replicates[[rep]]$seq.space.clusters$other$BIC))

    for(pc in PCs:1){
      coords[[pc]]<-cbind(coords[[pc]],replicates[[rep]]$seq.space.PCA$coordinates[,pc])
    }
    message(paste("bootstrap",rep))
  }
  message("combining bootstraps replicates")
  for(pc in PCs:1){
    coords.cor[[pc]]<-sqrt(cor(coords[[pc]],method=correlation)^2)
    diag(coords.cor[[pc]])<-NA
    coords.av.score[[pc]]<-mean(coords.cor[[pc]],na.rm = TRUE)
  }

  clust.opt <- which(colMeans(clust.fit)==max(colMeans(clust.fit)))
  score <- mean(coords.av.score)

  list(replicates = replicates,
       clust.fit  = clust.fit,
       clust.opt  = clust.opt,
       PC.scores  = coords.av.score,
       score      = score)
}

jackknife <- function(MSA,
                      res.prop   = res.prop2,
                      cys        = cys,
                      clusterPCs = 1:20,
                      clusters   = 1:15,
                      model      = "VVV",
                      sam        = 0.9,
                      reps       = 10,
                      PCs        = 3,
                      correlation= "pearson"){
  clust.fit     <- NULL
  replicates    <- NULL
  loads         <- NULL
  loads.cor     <- NULL
  loads.av.score<- NULL

  message(paste(reps, "jackknifed replicates using random",sam,"sequence samples"))
  for(rep in 1:reps){
    replicates[[rep]] <- PCA_MSA(MSA=MSA[sample(nrow(MSA), nrow(MSA)*sam),],
                                 res.prop   = res.prop,
                                 clusterPCs = clusterPCs,
                                 clusters   = clusters,
                                 model      = model,
                                 bootstrap  = NULL,
                                 jackknife  = NULL)

    clust.fit<-rbind(clust.fit,as.vector(replicates[[rep]]$seq.space.clusters$other$BIC))

    for(pc in PCs:1){
      loads[[pc]]<-cbind(loads[[pc]],replicates[[rep]]$seq.space.PCA$loadings[,pc])
    }
    message(paste("jackknife",rep))
  }
  message("combining jackknife replicates")
  for(pc in PCs:1){
    loads.cor[[pc]]<-sqrt(cor(loads[[pc]],method=correlation)^2)
    diag(loads.cor[[pc]])<-NA
    loads.av.score[[pc]]<-mean(loads.cor[[pc]],na.rm = TRUE)
  }

  clust.opt <- which(colMeans(clust.fit)==max(colMeans(clust.fit)))
  score<-mean(loads.av.score)

  list(replicates = replicates,
       clust.fit  = clust.fit,
       clust.opt  = clust.opt,
       PC.scores  = loads.av.score,
       score      = score)
}



subPCA_MSA <- function(SAPCA,
                       cluster       = 1,
                       subclusterPCs = 5,
                       subclusters   = 1:15,
                       model         = "VVV",
                       res.prop      = "inherit",
                       cys           = "inherit",
                       bootstrap     = NULL,
                       jackknife     = NULL,
                       sam           = 0.9,
                       repPCs        = 3,
                       correlation   = "pearson"){

  #################.
  # Define subset #
  #################.

  # Define subset for graph
  SUB = SAPCA$seq.space.clusters$classification==cluster
  #labels$M.AA.type.pc10.8Gs==1|2|3

  ###########################.
  # SUBSET PCA and clusters #
  ###########################.

  # Prepare data for PCA
  toPCA <- subset(SAPCA$numerical.alignment$MSA, subset = SUB)

  if(all(res.prop=="inherit")){
    res.prop <- SAPCA$numerical.alignment$res.prop
  }
  if(all(cys=="inherit")){
    cys <- SAPCA$call$cys
  }

  # PCA the subsetted data inclluding scaling
  subPCA <- PCA_MSA(MSA         = toPCA,
                    res.prop    = res.prop,
                    cys         = cys,
                    clusterPCs  = subclusterPCs,
                    clusters    = subclusters,
                    model       = model,
                    bootstrap   = bootstrap,
                    jackknife   = jackknife,
                    sam         = sam,
                    repPCs      = repPCs,
                    correlation = correlation)
  subPCA
}


# Analysis -------

topload <- function(SAPCA,
                    PC = 1,
                    n  = 20,
                    gapignore = FALSE){

  names     <- do.call(rbind, strsplit(gsub("\\.",":",rownames(SAPCA$seq.space.PCA$loadings)), ':'))
  consensus <- seqinr::consensus(SAPCA$numerical.alignment$MSA)
  combined  <- cbind(consensus,names,SAPCA$seq.space.PCA$loadings[,PC])
  sorted    <- combined[order(sqrt(SAPCA$seq.space.PCA$loadings[,PC]^2),decreasing = TRUE),]
  rownames(sorted) <- NULL
  colnames(sorted) <- c("consensus","MSA_col","property",paste("PC",PC,"_load",sep=""))
  if (gapignore){
    sorted <- sorted[sorted[,3]!="NOTGAP",]
  }
  if (n=="all") {
    sorted
  }else{
    head(sorted,n)
  }
}


loadingtable <- function(SAPCA,
                         PC        = 1,
                         seq       = seqinr::consensus(SAPCA$numerical.alignment$MSA),
                         magnitude = FALSE){

  input <- SAPCA$seq.space.PCA$loadings[,PC]

  if (magnitude==TRUE){
    input <- sqrt(input^2)
  }

  data <- matrix(data     = input,
                 nrow     = length(colnames(SAPCA$numerical.alignment$res.prop)),
                 byrow    = TRUE,
                 dimnames = list(colnames(SAPCA$numerical.alignment$res.prop),
                                   1:SAPCA$numerical.alignment$aln.len))
  sum             <- colSums(data)
  data2           <- rbind(data,sum)
  colnames(data2) <- paste(1:length(seq), "(", seq, ")", sep="")

  list(data       = data,
       annotated  = data2)
}


closest <- function (SAPCA,
                     sequence,
                     PC = 1:3,
                     n  = 10){

  coords     <- SAPCA$seq.space.PCA$coordinates
  centre     <- coords[sequence,PC]
  centre.m   <- matrix(rep(centre,nrow(coords)),
                       nrow  = nrow(coords),
                       byrow = TRUE)
  distances  <- SAPCA$seq.space.PCA$coordinates[,PC]-centre.m
  rootsquare <- sqrt(rowSums(distances^2))

  sorted     <- as.matrix(rootsquare[order(rootsquare)])
  colnames(sorted) <- "distance"

  head(sorted,n)
}


distance_matrix <- function (SAPCA,
                             PC     = 1:3,
                             subset = NULL){

  if (all(is.null(subset))){
    subset <- 1:nrow(SAPCA$seq.space.PCA$coordinates)
  }

  coords   <- SAPCA$seq.space.PCA$coordinates[subset,]
  output  <- NULL

  for (sequence in rownames(coords)) {
    centre     <- coords[sequence,PC]
    centre.m   <- matrix(rep(centre,nrow(coords)),
                         nrow  = nrow(coords),
                         byrow = TRUE)
    distances  <- coords[,PC]-centre.m
    rootsquare <- as.matrix(sqrt(rowSums(distances^2)))
    colnames(rootsquare) <- sequence
    output    <- cbind(output,rootsquare)
  }

  output
}



PyMol_columns <- function (SAPCA,
                           template,
                           PC=1,
                           write=FALSE){
  # Which columns of the template are not gaps
  structural <- SAPCA$numerical.alignment$MSA[template,]!="-"

  # For the temaplte columns, tabulate the residue loadings for the elected PC
  output <- 0*loadingtable(SAPCA)$data[,structural]
  for(i in PC){
    output <- output + loadingtable(SAPCA,PC = i,magnitude = T)$data[,structural]
  }
  output <- rbind(10*output,10*colSums(output),SAPCA$numerical.alignment$MSA[template,structural])
  if (write!=FALSE){
    filename <- paste(name,".Structureloadings",".PC",PC,".(",template,").csv", sep="")
    if (write!=TRUE){filename <- write}
    write.csv(t(output),filename)
  }
  output
}


# Adding sequences -------

seq.MSA.add <- function(SAPCA,
                        sequence,
                        SAPCAname=NULL){
  sequence <- casefold(sequence,upper=TRUE)
  MSA   <- SAPCA$numerical.alignment$MSA
  MSA2  <- as.AAstringSet(MSA,degap = TRUE)
  seqs  <- nrow(MSA)
  seq   <- as.AAstring(sequence, degap=FALSE)
  seq.d <- as.AAstring(sequence, degap=TRUE)
  BLOSUM40 <- read.blosum()

  aln.all <- Biostrings::pairwiseAlignment(MSA2,
                                           seq.d,
                                           substitutionMatrix = BLOSUM40,
                                           gapOpening   = 0,
                                           gapExtension = 8,
                                           scoreOnly    = TRUE)

  # Max possible similarity score
  aln.limit <- Biostrings::pairwiseAlignment(seq.d,
                                             seq.d,
                                             substitutionMatrix = BLOSUM40,
                                             gapOpening   = 0,
                                             gapExtension = 8,
                                             scoreOnly    = TRUE)

  # Similarity score as percentage of max
  aln.hit.score <- max(aln.all)/aln.limit

  # The sequence of the best matching member of the database
  aln.hit.num  <- which(aln.all==max(aln.all))[1]
  aln.hit.name <- SAPCA$numerical.alignment$seq.names[aln.hit.num]
  aln.hit.seq  <- paste(as.AAstring(MSA[aln.hit.num,],degap = 1))

  # Use "*" to indicate gaps in the best reference sequence (aln.hit)
  aln.hit <- gsub("-","O",as.AAstring(MSA[aln.hit.num,]))
  tomsa   <- as.AAstringSet(rbind(as.character(aln.hit),
                                  as.character(seq.d)),
                            degap=TRUE)
  aln.add <- msa(tomsa,
                 substitutionMatrix = BLOSUM40,
                 gapOpening   = 0,
                 gapExtension = 1)

  # Has the new sequence introduced exrta gaps into the hit sequence alignement?
  aln.hit.orig         <- as.AAstring(MSA[aln.hit.num,])
  aln.hit.new          <- as.AAstring(as.character(aln.add)[1])
  aln.hit.seq.aln.orig <- unlist(strsplit(as.character(aln.hit.orig),""))
  aln.hit.seq.aln.new  <- unlist(strsplit(as.character(aln.hit.new),""))

  if(as.character(aln.hit.orig)!=as.character(aln.hit.new)){
    message(paste(sum(aln.hit.seq.aln.new=="-"),
                "residues ofthe new sequence were not alignable to the",
                SAPCAname,
                "reference MSA so have been ignored"))
  }

  # Alignment addition as matrix
  aln.add.mat <- unlist(as.matrix(aln.add),"")

  # Unmatchable resiues removed from aligned sequence
  aln.add2 <- aln.add.mat[2,][aln.add.mat[1,]!="-"]
  aln.add3 <- paste(as.AAstring(aln.add2))

  # Gaps in the hit sequence (original and newly aligned)
  gaps.orig       <- unlist(strsplit(as.character(aln.hit.orig),"[A-Z]"))
  gaps.count.orig <- nchar(gaps.orig)
  gaps.lead.orig  <- gaps.count.orig[1]
  gaps.trail.orig <- gaps.count.orig[length(gaps.count.orig)]

  gaps.new        <- unlist(strsplit(as.character(gsub("-","",aln.hit.new)),"[^O]"))
  gaps.count.new  <- nchar(gaps.new)
  gaps.lead.new   <- gaps.count.new[1]
  gaps.trail.new  <- gaps.count.new[length(gaps.count.new)]

  if(length(gaps.count.orig)>length(gaps.count.new)){
    gaps.count.new <- append(gaps.count.new,0)
  }

  gaps.discrep     <- suppressWarnings(rbind(gaps.count.orig,gaps.count.new))
  gaps.discrep.num <- gaps.discrep[1,]-gaps.discrep[2,]
  gaps.lead.add    <- gaps.discrep.num[1]
  gaps.trail.add   <- gaps.discrep.num[length(gaps.discrep.num)]

  # New alignment
  aln.add4  <- c(rep("-",gaps.lead.add),
                 aln.add2,
                 rep("-",gaps.trail.add))

  original      <- as.character(seq.d)
  alignable     <- aln.add3
  
  tomsa2 <- as.AAstringSet(rbind(original,
                                 alignable),
                           degap=TRUE)
 
   seq.alignable <- msa::msa(tomsa2,
                             substitutionMatrix = BLOSUM40,
                             gapOpening   = 0,
                             gapExtension = 1)
  length(aln.add4)==ncol(MSA)

  query     <- aln.add4

  aln.final <- rbind(query,MSA)
  output    <- list(MSA             = aln.final,
                    aln.hit.name    = aln.hit.name,
                    aln.hit.seq     = aln.hit.seq,
                    aln.hit.score   = aln.hit.score,
                    aln.all.score   = aln.all/aln.limit,
                    seq.alignable   = seq.alignable)
  output
}



seq.rotate <- function(SAPCA,newseq){

  res.props  <- colnames(SAPCA$numerical.alignment$res.prop)
  prop.means <- SAPCA$numerical.alignment$prop.means
  prop.vars  <- SAPCA$numerical.alignment$prop.vars
  # Align new sequence with MSA

  # Numericise new sequence
  newseq.num <- numericise_MSA(MSA      = newseq$MSA[c("query",newseq$aln.hit.name),],
                               res.prop = SAPCA$numerical.alignment$res.prop)

  # Scale new sequnce using same scaling as SAPCA (gaps as "NA")
  newseq.scale.stack <- NULL
  for (x in 1:length(res.props)) {
    newseq.scale.stack[[res.props[x]]] <- (newseq.num$MSA.num.stack[[res.props[x]]]- prop.means[x]) /
      sqrt(prop.vars[x])
  }

  # Reflow into single wide matrix
  newseq.scale.wide <- newseq.scale.stack[[1]]
  for (x in 2:length(res.props)) {
    newseq.scale.wide <- cbind(newseq.scale.wide, newseq.scale.stack[[res.props[x]]])
  }

  # Replace gaps (currently "NA") with column average of the scaled SAPCA
  # Create na.colmean function
  gapvalues  <- colMeans(SAPCA$numerical.alignment$MSA.scale.wide)
  na.replace <- function(x,y){
    x[is.na(x)] <- y
    x
  }
  newseq.scale.wide.g <- NULL
  for(i in 1:ncol(newseq.scale.wide)){
    newseq.scale.wide.g <- cbind(newseq.scale.wide.g,
                                 na.replace(newseq.scale.wide[,i],gapvalues[i]))
  }

  # Rotate scaled sequence into same space as SAPCA
  newseq.rot <- scale(newseq.scale.wide.g,
                      SAPCA$seq.space.PCA$centre,
                      SAPCA$seq.space.PCA$scale) %*% SAPCA$seq.space.PCA$loadings

  # Output
  output <- list(seq       = newseq,
                 seq.num   = newseq.num$MSA.num.wide,
                 seq.scale = newseq.scale.wide.g,
                 seq.rot   = newseq.rot)
  output
}



seq.clust.add <- function(SAPCA,newseq.r){

  SAPCA.c  <- mclustrev(SAPCA)
  newseq.c <- mclust::predict.Mclust(SAPCA.c,newseq.r$seq.rot[,SAPCA$call$clusterPCs])
  newseq.c
}



mclustrev <- function(SAPCA){
  output <- SAPCA$seq.space.clusters$other

  output$classification <- SAPCA$seq.space.clusters$classification
  output$G              <- SAPCA$seq.space.clusters$optimal
  output$z              <- SAPCA$seq.space.clusters$likelihoods

  class(output) <- "Mclust"
  output
}



seq.SAPCA.add <- function (SAPCA,newseq.r,newseq.c){
  output <- SAPCA
  output$numerical.alignment$seq.names      <- rownames(newseq.r$seq$MSA)
  output$numerical.alignment$MSA            <- newseq.r$seq$MSA
  output$numerical.alignment$MSA.num.wide   <- rbind(newseq.r$seq.num[1,],
                                                     SAPCA$numerical.alignment$MSA.num.wide)
  output$numerical.alignment$MSA.scale.wide <- rbind(newseq.r$seq.scale[1,],
                                                     SAPCA$numerical.alignment$MSA.scale.wide)
  output$numerical.alignment$MSA.num.stack  <- NULL
  output$numerical.alignment$MSA.scale.stack<- NULL

  output$seq.space.PCA$coordinates          <- rbind(newseq.r$seq.rot[1,],
                                                     SAPCA$seq.space.PCA$coordinates)
  output$seq.space.clusters$likelihoods     <- rbind(newseq.c$z[1,],
                                                     SAPCA$seq.space.clusters$likelihoods)
  output$seq.space.clusters$classification  <- c(newseq.c$classification[1],
                                                 SAPCA$seq.space.clusters$classification)

  rownames(output$numerical.alignment$MSA.num.wide)[1]   <- "query"
  rownames(output$numerical.alignment$MSA.scale.wide)[1] <- "query"
  rownames(output$seq.space.PCA$coordinates)[1]          <- "query"
  rownames(output$seq.space.clusters$likelihoods)[1]     <- "query"

  output
}



seq.add.full <- function (SAPCA,sequence,SAPCAname=NULL){
  sequence <- casefold(sequence,upper=TRUE)
  newseq   <- seq.MSA.add(SAPCA,sequence,SAPCAname)
  newseq.r <- seq.rotate(SAPCA,newseq)
  newseq.c <- seq.clust.add(SAPCA,newseq.r)
  SAPCA2   <- seq.SAPCA.add(SAPCA,newseq.r,newseq.c)

  SAPCA2
}


# Input/output --------

###############.
# Prepare MSA #
###############.

read.MSA <- function(MSA){
  # Load sequence MSA
  # if a matrix, can be used straight away
  # if raw fasta file, use seqinr to convert to data frame
  if (!is.matrix(MSA)){
    MSA       <- data.frame(seqinr::read.fasta(MSA,set.attributes=FALSE))
  }
  # if a data frame, convert to matrix
  if (is.data.frame(MSA)){
    MSA       <- as.matrix(t(toupper(as.matrix(MSA))))
  }

  MSA
}

##############################.
# Prepare residue properties #
##############################.

read.res.prop <- function(res.prop, cys=TRUE){
  # Load residue properties
  # if a data frame, can be used straight away
  # if a matrix, convert to data frame
  if (is.matrix(res.prop)){
    res.prop <- data.frame(res.prop)
  }
  # if a csv file, convert to data frame
  if (!is.data.frame(res.prop)){
    res.prop <- read.csv(res.prop)
  }
  if (all(row.names(res.prop)==1:nrow(res.prop))){
    # Format rownames and order
    row.names(res.prop) <- c(as.matrix(res.prop[,1]))           # first column to row names
    res.prop            <- res.prop[,-1]                        # remove original first column
    res.prop            <- res.prop[,order(colnames(res.prop))] # order columns alphabetically
    # Append X and gap rows (if X not already defined)
    if (length(grep("X", rownames(res.prop)))==0){
      res.prop          <- rbind(res.prop, colMeans(res.prop))
      rownames(res.prop)[nrow(res.prop)] <- "X"
    }
    res.prop            <- rbind(res.prop, rep(NA,ncol(res.prop)))
    rownames(res.prop)[nrow(res.prop)]   <- "-"
    # Add cystene columns (gaps = 0)
    if (cys==TRUE){
      res.prop[ncol(res.prop)+1]                    <- 0
      colnames(res.prop)[ncol(res.prop)]            <- "CYS"
      res.prop[grep("C", rownames(res.prop)),"CYS"] <- 1
    }
    # Append presence/absence columns
    res.prop[ncol(res.prop)+1]         <- c(rep(1,nrow(res.prop)-1),0)
    colnames(res.prop)[ncol(res.prop)] <- "NOTGAP"
  }
  res.prop
}



read.blosum <- function(file="C:\\Users\\T\\OneDrive\\0-Sequences\\2-PCA\\0-Raw data and scalers\\0 - BLOSUM40.csv"){
  BLOSUM40 <- read.csv(file)
  BLOSUM40.names <- BLOSUM40[,1]
  BLOSUM40 <- BLOSUM40[,-1]
  BLOSUM40 <- as.matrix(BLOSUM40)
  rownames(BLOSUM40)<-BLOSUM40.names
  colnames(BLOSUM40)<-BLOSUM40.names
  BLOSUM40
}



as.fasta <- function(matrix,degap=FALSE,decolgap=FALSE,write=FALSE,name=NULL){
  
  # Remove empty columns
  if(decolgap){
    matrix<-matrix[,colMeans(matrix=="-")!=1]
  }
  
  if(is.list(matrix)){
    # Convert alignment matrix to list of strings
    seqs <- unlist(as.vector(lapply(matrix,paste, collapse = "")))
    names <- paste0(">",names(matrix))
  }else{
    # Convert alignment matrix to list of strings
    seqs  <- do.call("paste",c(data.frame(matrix),sep=""))
    names <- paste(">",row.names(matrix),sep="")
    
    # If just one sequence, this is how to name it
    if(is.null(dim(matrix))){
      names <- ">sequence"
      if(!is.null(name)){
        names <- paste0(">",name)
      }
      seqs  <- paste(matrix,collapse="")
    }
  }
  
  # Degap sequences
  if (degap){
    seqs <- gsub("-","",seqs)
  }
  
  # Interleave names and sequences
  ord1 <- 2*(1:length(names))-1
  ord2 <- 2*(1:length(seqs))
  
  # Output
  if (write==FALSE){
    cat(c(names,seqs)[order(c(ord1,ord2))], sep = "\n")
  }else{
    if (!grepl(".fa",write,ignore.case=TRUE)){
      write<-paste(write,".fa",sep="")
    }
    cat(c(names,seqs)[order(c(ord1,ord2))], sep = "\n", file = write)
  }
}


as.AAstring<-function(string, degap=FALSE){
  string <- paste(string,collapse="")
  if(degap==TRUE){
    string<-gsub("-","",string)
  }
  output <- Biostrings::AAString(string)
  output
}



as.AAstringSet<-function(MSA, degap=FALSE){
  MSA <- apply(MSA,1,paste,collapse="")
  if(degap==TRUE){
    MSA<-gsub("-","",MSA)
  }
  output <- Biostrings::AAStringSet(MSA)
  output
}



percent <- function(x, digits = 1, format = "f", ...) {
  paste0(formatC(100 * x, format = format, digits = digits, ...), "%")
}


