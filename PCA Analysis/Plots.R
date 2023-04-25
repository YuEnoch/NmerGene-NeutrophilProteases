#########.
# PLOTS #
#########.

# Basic plots ---------

plot_modelfit <- function(SAPCA,top=0.001,legend=TRUE){
  # Model plots     
  if(!is.null(SAPCA$bootstrap.reps)){
    data    <- as.matrix(colMeans(SAPCA$bootstrap.reps$clust.fit))
    data.se <- apply(SAPCA$bootstrap.reps$clust.fit, 2, sd, na.rm = 1)/
               sqrt(nrow(SAPCA$bootstrap.reps$clust.fit))
    opt     <- SAPCA$bootstrap.reps$clust.opt
  }
  if(!is.null(SAPCA$jackknife.reps)){
    data    <- as.matrix(colMeans(SAPCA$jackknife.reps$clust.fit))
    data.se <- apply(SAPCA$jackknife.reps$clust.fit, 2, sd, na.rm = 1)/
               sqrt(nrow(SAPCA$jackknife.reps$clust.fit))
    opt     <- SAPCA$jackknife.reps$clust.opt
  }
  if(is.null(SAPCA$bootstrap.reps) && is.null(SAPCA$jackknife.reps)){
    data    <- as.matrix(SAPCA$seq.space.clusters$other$BIC)
    data.se <- NULL
    opt     <- SAPCA$seq.space.clusters$optimal
  }
  
  optname <- SAPCA$seq.space.clusters$other$modelName
  optpos  <- which(data==max(data,na.rm = TRUE))
  opts    <- (data > quantile(data,1-top,na.rm=TRUE))*row(data)
  opts[opts==0] <- NA
  labels  <- rep(NA,length(data))
  if(ncol(data)==1){
    labels[optpos] <- opt
  }else{
    labels[optpos] <- paste(opt,optname)
    matplot(data[,],
            type = "b",
            pch  = 1,
            ylim = c(min(data,na.rm = TRUE),max(data,na.rm = TRUE)*0.9),
            xlim = c(1,nrow(data)+1),
            xlab = "Number of groups",
            ylab = "Goodness of fit")
  }
  
    plot(data[,],
          type = "l",
          pch  = 1,
          ylim = c(min(data,na.rm = TRUE),max(data,na.rm = TRUE)*0.9),
          xlim = c(1,nrow(data)+1),
          xlab = "Number of groups",
          ylab = "Goodness of fit")
  
  if(top>0.001){
    lines(x=0:(nrow(data)+2),
          y=rep(quantile(data,1-top,na.rm=TRUE,type=2),nrow(data)+3),
          lty=3,col="red")
  }else{
    lines(x=rep(opt,2),
          y=c(min(data,na.rm = TRUE)*1.5,max(data,na.rm = TRUE)),
          lty=3,col="red")
  }
  text(x=opts*(!is.na(opts)),y=data*(!is.na(opts)),labels = labels,pos=3)
  points(x=opts[!is.na(opts)],y=data[!is.na(opts)], col = "red",pch=16)
  if(all(ncol(data)>1,legend==TRUE)){
    legend(x      = "bottomright",
           legend = colnames(data),
           col    = 1:ncol(data),
           ncol   = round(ncol(data)/4+0.5),
           pch    = 1,
           bty    = "n")
  }
  if(!is.null(data.se)){
    mclust::errorBars(1:length(data),
                      upper = data + data.se,
                      lower = data - data.se,
                      width = 0.04)
  }
}



plot_scree <- function(SAPCA){
  # Scree plot of component significance
  barplot(SAPCA$seq.space.PCA$stdev[1:15],               # first 15 principal components
          xlab = "Principal component",        # x label
          ylab = "Variance",                  # y label
          main = "Principal components")      # title
}



plot_heatmap <- function(SAPCA,
                         PC = 1:3){
  
  distances <- distance_matrix(SAPCA,PC=PC)
  heatmap(distances)
}  



plot_network <- function(SAPCA,
                         PC      = 1:3,
                         toplink = 20,
                         col     = "cluster",
                         size    = 4){
  
  if(all(col=="cluster")){
    colour <- colours[SAPCA$seq.space.clusters$classification]
  }else{
    colour <- col
  }

  distances <- distance_matrix(SAPCA, PC=PC)
  distances[distances>quantile(distances, prob=toplink/100)] <- 0
  g <- igraph::graph.adjacency(distances, mode="undirected", diag=FALSE, weighted=TRUE)
  plot(g,
       #edge.width   = linewidth*max(E(g)$weight)/E(g)$weight,
       edge.curved  = 0.2,
       vertex.size  = size,
       vertex.label = NA,
       vertex.color = colour)
}   



plot_loadings_matrix <- function(SAPCA,
                                 PC        = 1,
                                 topload   = 100,
                                 magnitude = FALSE,
                                 occupancy = FALSE){
  
  d  <- loadingtable(SAPCA,PC=PC,magnitude=magnitude)$data
  d[sqrt(d^2) < quantile(sqrt(d^2),prob=1-topload/100)] <- 0
  
  if(all(occupancy==FALSE)){
    occupancy <- NULL
    div.line  <- data.frame(x=0,y=0,xend=0,yend=0)
  }else{
    occupancy <- colMeans(SAPCA$numerical.alignment$MSA.num.stack$NOTGAP)*max(d)*0.8
    div.line  <- data.frame(x    = 0.5,
                            y    = nrow(d)+0.5,
                            xend = ncol(d)+0.5,
                            yend = nrow(d)+0.5)
  }
  
  d  <- rbind(d,occupancy)
  md <- reshape2::melt(d)
  colnames(md) <- c('Property',"Residue_number",'Loading')
  
  p <-  ggplot2::ggplot(md)                                                                       +
        ggplot2::geom_tile(ggplot2::aes_string(x="Residue_number", y="Property", fill="Loading")) +
        ggplot2::scale_x_continuous(expand = c(0, 0))                                             +
        ggplot2::scale_y_discrete  (expand = c(0, 0))                                             +
        ggplot2::theme_classic()                                                                  +
        ggplot2::theme(panel.border = ggplot2::element_rect(colour = "black", fill=NA, size=1))   +
        ggplot2::geom_segment(data=div.line, ggplot2::aes(x,y,xend=xend, yend=yend), size=0)      +
        ggplot2::scale_fill_gradient2(low      = "darkred",
                                      high     = "darkblue",
                                      mid      = "white",
                                      na.value = "grey50")
  p
}



plot_loadings_categories <- function(SAPCA,
                                     PC      = 1,
                                     topload = 10){

  md <- reshape2::melt(loadingtable(SAPCA,PC,magnitude=TRUE)$data)
  colnames(md) <- c('Property',"Residue_number",'Loading')
  
  md[md$Loading < quantile(md$Loading,prob=1-topload/100),"Loading"] <- 0
  md[md$Loading>0 ,"Loading"] <- as.character(md[md$Loading>0 ,"Property"])
  
  p <-  ggplot2::ggplot(md)                                                                       +
        ggplot2::geom_tile(ggplot2::aes_string(x="Residue_number", y="Property", fill="Loading")) +
        ggplot2::scale_x_continuous(expand = c(0, 0))                                             +
        ggplot2::scale_y_discrete  (expand = c(0, 0))                                             +
        ggplot2::theme_classic()                                                                  +
        ggplot2::theme(panel.border = ggplot2::element_rect(colour = "black", fill=NA, size=1))   +
        ggplot2::scale_fill_manual(guide=FALSE,
                                   values=c("white",
                                            "red",
                                            "blue",
                                            "darkgreen",
                                            "purple",
                                            "grey40"))   
  p
}


# 3D plots ---------

plot_3Dclusters <- function(SAPCA,
                            plotPCs    = 1:3,
                            col        = "cluster",
                            radius     = 1,
                            labels     = NULL,
                            labeltext  = NULL,
                            write      = FALSE,
                            axeslabels = "default"){
  if(!is.null(SAPCA$seq.space.PCA$coordinates)){
    data <- SAPCA$seq.space.PCA$coordinates
  }else{
    data <- SAPCA    
  }
  
  if(all(col=="cluster")){
    colour <- SAPCA$seq.space.clusters$classification
  }else{
    colour <- col
  }
  # Calculate radius size
  rad <- (range(SAPCA$seq.space.PCA$coordinates[,plotPCs])[2] -
          range(SAPCA$seq.space.PCA$coordinates[,plotPCs])[1]) /
          120
  rad <- rad*radius
  
  if(all(axeslabels=="default")){
    axes <- paste("PC",plotPCs,sep="")
  }else{
    axes <- axeslabels
  }
  if(is.null(axeslabels)){
    axes <- c("","","")
  }
    
  # Plot model-based clusters in 3D
  rgl::plot3d(data[,plotPCs],
              col      = colour,      # colour by clusters
              specular = "black",     # matte lighting
              type     = "s",         # "p" is points, "s" is spheres
              radius   = rad,         # sphere radius if using spheres
              size     = 4,           # point size
              axes     = FALSE,       # draw axes separately
              xlab     = axes[1],
              ylab     = axes[2],
              zlab     = axes[3])       
  # Draw axes
  if(write!=FALSE){
    rgl::axes3d(color = "black", labels = FALSE)                       
  }else{
    rgl::axes3d(color = "black", alpha=0.5, labels = FALSE) 
  }
  
  for (NAME in labels){
    # Which point will be labelled
    if(is.numeric(labels) | is.logical(labels)){
      SUB <- NAME    # Label based on its order number
    }else{
      SUB <- row.names(SAPCA$seq.space.PCA$coordinates)==NAME  # Label based on its row.name
    }
    
    # What is the label text
    if(is.null(labeltext)){
      TEXT <- NAME
      if(is.numeric(labels) | is.logical(labels)){
        TEXT <- SAPCA$numerical.alignment$seq.names[NAME]
      }
    }else{
      TEXT <- labeltext[labels==NAME] 
    }
    
    rgl::text3d(SAPCA$seq.space.PCA$coordinates[SUB,plotPCs], 
                text      = paste('---',TEXT),   # data label text
                font      = 2,                   # bold
                color     = "black",             # colour
                adj       = -rad/2)              # offset
  }
  
  # Write html for interactive data
  if(write!=FALSE){
    rgl::writeWebGL(write)                          
  }
}



plot_3Dbiplot <- function(SAPCA,
                          plotPCs   = 1:3,
                          col       = "property",
                          radius    = 1,
                          labels    = NULL,
                          labeltext = NULL,
                          write     = FALSE){
  if(all(col=="property")){
    colour <- NULL
    for (x in 1:ncol(SAPCA$numerical.alignment$res.prop)){
      colour <- append(colour,rep(x,SAPCA$numerical.alignment$aln.len))
    }
  }else{
    colour <- col
  }
  # Calculate radius size
  rad <- (range(SAPCA$seq.space.PCA$loadings[,1])[2] - 
          range(SAPCA$seq.space.PCA$loadings[,1])[1]) /
          100
  rad <- rad*radius
  
  # Plot model-based clusters in 3D
  rgl::plot3d(SAPCA$seq.space.PCA$loadings[,plotPCs],
              col      = colour,      # colour by clusters
              specular = "black",     # matte lighting
              type     = "s",         # "p" is points, "s" is spheres
              radius   = rad,         # sphere radius if using spheres
              size     = 4,           # point size
              axes     = FALSE)       # draw axes separately
  # Draw axes
  if(write!=FALSE){
    rgl::axes3d(color = "black", labels = FALSE)                       
  }else{
    rgl::axes3d(color = "black", alpha=0.5, labels = FALSE) 
  }
  
  for (NAME in labels){
    
    # Which point will be labelled
    if(is.numeric(labels) | is.logical(labels)){
      SUB <- NAME    # Label based on its order number
    }else{
      SUB <- row.names(SAPCA$seq.space.PCA$loadings)==NAME  # Label based on its row.name
    }
    
    # What is the label text
    if(is.null(labeltext)){
      TEXT <- NAME
      if(is.numeric(labels) | is.logical(labels)){
        TEXT <- rownames(SAPCA$seq.space.PCA$loadings)[NAME]
      }
    }else{
      TEXT <- labeltext[labels==NAME]
    }
    
    rgl::text3d(SAPCA$seq.space.PCA$loadings[SUB,plotPCs], 
              text      = paste('---',TEXT),   # data label text
              font      = 2,                   # bold
              color     = "black",             # colour
              adj       = -rad/2)              # offset
  }
  
  # Write html for interactive data
  if(write!=FALSE){
    rglwidget::.writeWebGL(write)                          
  }
}



plot_3Dclosest <- function(SAPCA,
                           sequence,
                           plotPCs    = 1:3,
                           measurePCs = 1:3,
                           radius     = 1,
                           n          = 10,
                           write      = FALSE,
                           axeslabels = "default"){

  temp <- SAPCA$numerical.alignment$MSA[,1]
  temp[1:length(temp)]                    <- "white"
  temp[rownames(closest(SAPCA,sequence,PC = measurePCs, n = n))] <- "red"
  temp[sequence]                          <- "black"

  plot_3Dclusters(SAPCA,
                  radius     = radius,
                  plotPCs    = plotPCs,
                  labels     = sequence,
                  col        = temp,
                  write      = write,
                  axeslabels = axeslabels)
}



plot_3Drotation <- function(write="SeqSpace", axis = c(0,0,1), fps=24, duration=8, loop=TRUE, clean=TRUE){
  
  folder <- paste(write, "3Dplot.rotation", sep=".")
  dir.create(folder)
  
  rgl::movie3d(rgl::spin3d(axis = axis,
                           rpm  = 25/duration),
               duration = fps*duration/10,
               type     = "gif",
               clean    = FALSE,
               convert  = FALSE,
               top      = FALSE,
               dir      = paste(getwd(),folder,sep="/"))
  
  frames     <- duration*fps
  padding    <- if(stringr::str_length(frames)<=3){3}else{stringr::str_length(frames)}
  framenames <- NULL
  for(frame in 0:frames){
    framenames <- append(framenames,paste(folder,
                                          "/movie",
                                          stringr::str_pad(frame,
                                                           padding,
                                                           pad = "0"),
                                          ".png",sep=""))
  }
  
  # Animation settings
  animation::ani.options(loop = loop)      # Looping
  animation::ani.options(interval = 1/fps) # Framerate
  
  #Make 100-frame sections
  sections <- round(frames/100+0.5)
  sectionnames <- paste(folder,"/animation_section",1:sections,".gif",sep="")
  for(section in 1:sections){
    animation::im.convert(files  = na.exclude(framenames[((section-1)*100)+(1:100)]),
                          output = sectionnames[section])
  }
  # Combine sections to final video
  animation::im.convert(files  = sectionnames,
                        output = paste(folder,".gif",sep=""))
  
  # Clean up
  unlink(sectionnames)
  if(clean==TRUE){
    unlink(folder, recursive = TRUE)
  }
}



plot_overlay_3Dnetwork <- function(SAPCA,
                                   toplink = 5,
                                   plotPCs = 1:3,
                                   netPCs  = "default",
                                   subset  = NULL,
                                   col     = "black",
                                   alpha   = 0.2,
                                   lwd     = 1){

  if(all(netPCs=="default")){
    netPCs <- plotPCs
  }
  
  if(all(is.null(subset))){
    subset <- 1:nrow(SAPCA$seq.space.PCA$coordinates)
  }
  
  distances  <- distance_matrix(SAPCA, PC=netPCs, subset=subset)
  mat        <- distances<=quantile(distances, prob=toplink/100)
  list       <- reshape2::melt(mat)[as.vector(mat),]
  coords     <- SAPCA$seq.space.PCA$coordinates[subset,plotPCs]
  line.start <- coords[list[,1],]
  line.end   <- coords[list[,2],]
  
  lines.to.add <- matrix(data     = rbind(t(line.start),
                                          t(line.end)),
                         ncol     = 3,
                         byrow    = TRUE,
                         dimnames = list(c(rbind(rownames(line.start),
                                                 rownames(line.end))),
                                         c("x","y","z")))
  
  rgl::segments3d(lines.to.add,
                  col   = col,
                  alpha = alpha,
                  lwd   = lwd)
}



plot_overlay_3Dtree <- function (SAPCA     = SAPCA,
                               tree      = "nj",
                               plotPCs   = 1:3,
                               njPCs     = 1:3,
                               ancestors = NULL,
                               col       = "black",
                               alpha     = 0.2,
                               lwd       = 1){

  #zero-length distances have to be replaces with some positive integer for FastAnc
  if(all(tree=="nj")){
    dist          <- distance_matrix(SAPCA,PC = njPCs)
    dist[dist==0] <- min(dist[dist!=0])/10
    tree          <- ape::nj(dist)
  }
  if(is.null(attributes(tree)[[2]])){
    tree <- ape::read.tree(tree)
  }
  
  names <- tree$tip.label
  
  #format coordinates for tips
  X <- SAPCA$seq.space.PCA$coordinates[names,plotPCs]
  
  #if tree  has no edge lengths, then 
  if(is.null(tree$edge.length)){
    tree$edge.length <- rep(1,nrow(tree$edge))
  }
  #remove 0-length edges
  tree$edge.length[tree$edge.length<=0] <- min(tree$edge.length[tree$edge.length>0])
  
  
  #ancestral node 3d locations
  if(is.null(ancestors)){ 
    A <- apply(X,
               2,
               function(x, tree) phytools::fastAnc(tree, x), 
               tree = tree)
  }else{
    A <- ancestors[as.character(1:tree$Nnode + length(tree$tip.label)),]
  }
  
  #prepare data order
  x <- y <- z <- matrix(NA, nrow(tree$edge), 2)
  X <- X[tree$tip.label, ]
  #format tip locations
  for (i in 1:length(tree$tip.label)) {
    x[tree$edge[, 2] == i, 2] <- X[i, 1]
    y[tree$edge[, 2] == i, 2] <- X[i, 2]
    z[tree$edge[, 2] == i, 2] <- X[i, 3]
  }
  #format node locations
  for (i in length(tree$tip.label) + 1:tree$Nnode) {
    x[tree$edge[, 1] == i, 1] <- x[tree$edge[, 2] == i, 2] <- A[as.character(i),1]
    y[tree$edge[, 1] == i, 1] <- y[tree$edge[, 2] == i, 2] <- A[as.character(i),2]
    z[tree$edge[, 1] == i, 1] <- z[tree$edge[, 2] == i, 2] <- A[as.character(i),3]
  }
  
  # #line colours
  # tip.colour  <- SAPCA$seq.space.clusters$classification[names]
  # anc.colour  <- round(phytools::fastAnc(tree2, C),0)
  # all.colour  <- c(tip.colour,anc.colour)
  # outer       <- tree$edge[,2]
  # outer[tree$edge[,1]<=tree$edge[,2]] <- tree$edge[tree$edge[,1]<=tree$edge[,2],1]
  # line.colour <- all.colour[outer]
  
  #format line locations
  line.start   <- cbind(x[,1],y[,1],z[,1])
  line.end     <- cbind(x[,2],y[,2],z[,2])
  lines.to.add <- matrix(data     = rbind(t(line.start),
                                          t(line.end)),
                         ncol     = 3,
                         byrow    = TRUE,
                         dimnames = list(c(rbind(rownames(line.start),
                                                 rownames(line.end))),
                                         c("x","y","z")))
  
  #plot phylogeny lines
  rgl::segments3d(lines.to.add,
                  col   = col,
                  alpha = alpha,
                  lwd   = lwd)
}


plot_overlay_3Dlabel <- function(SAPCA,
                                 plotPCs = 1:3){
  selected     <- rgl::select3d()
  selected.set <- selected(SAPCA$seq.space.PCA$coordinates[,plotPCs])
  if(sum(selected.set)!=0){
    as.fasta(SAPCA$numerical.alignment$MSA[selected.set,],
             decolgap=TRUE)
    
    
    for (NAME in SAPCA$numerical.alignment$seq.names[selected.set]){
      # Which point will be labelled
      SUB <- NAME   
      # What is the label text
      TEXT <- NAME
      
      rad <- (range(SAPCA$seq.space.PCA$coordinates[,plotPCs])[2] -
              range(SAPCA$seq.space.PCA$coordinates[,plotPCs])[1]) /
              100
      
      rgl::text3d(SAPCA$seq.space.PCA$coordinates[SUB,plotPCs], 
                  text      = paste('---',TEXT),   # data label text
                  font      = 2,                   # bold
                  color     = "black",             # colour
                  adj       = -rad/2)              # offset
    }
  }
}
  
