library("phyloseq")
library("data.table")
library("ggplot2")
library(maps)
library(mapdata)
library(scales)
library(mapproj)
library(intrval)

#from here https://github.com/joey711/phyloseq/issues/418
fast_melt = function(physeq){
  # supports "naked" otu_table as `physeq` input.
  otutab = as(otu_table(physeq), "matrix")
  if(!taxa_are_rows(physeq)){otutab <- t(otutab)}
  otudt = data.table(otutab, keep.rownames = TRUE)
  setnames(otudt, "rn", "taxaID")
  # Enforce character taxaID key
  otudt[, taxaIDchar := as.character(taxaID)]
  otudt[, taxaID := NULL]
  setnames(otudt, "taxaIDchar", "taxaID")
  # Melt count table
  mdt = melt.data.table(otudt, 
                        id.vars = "taxaID",
                        variable.name = "SampleID",
                        value.name = "count")
  # Remove zeroes, NAs
  #mdt <- mdt[count > 0][!is.na(count)]
  # Calculate relative abundance
  mdt[, RelativeAbundance := count / sum(count), by = SampleID]
  if(!is.null(tax_table(physeq, errorIfNULL = FALSE))){
    # If there is a tax_table, join with it. Otherwise, skip this join.
    taxdt = data.table(as(tax_table(physeq, errorIfNULL = TRUE), "matrix"), keep.rownames = TRUE)
    setnames(taxdt, "rn", "taxaID")
    # Enforce character taxaID key
    taxdt[, taxaIDchar := as.character(taxaID)]
    taxdt[, taxaID := NULL]
    setnames(taxdt, "taxaIDchar", "taxaID")
    # Join with tax table
    setkey(taxdt, "taxaID")
    setkey(mdt, "taxaID")
    mdt <- taxdt[mdt]
  }
  return(mdt)
}




summarize_taxa = function(physeq, Rank, GroupBy = NULL){
  Rank <- Rank[1]
  if(!Rank %in% rank_names(physeq)){
    message("The argument to `Rank` was:\n", Rank,
            "\nBut it was not found among taxonomic ranks:\n",
            paste0(rank_names(physeq), collapse = ", "), "\n",
            "Please check the list shown above and try again.")
  }
  if(!is.null(GroupBy)){
    GroupBy <- GroupBy[1]
    if(!GroupBy %in% sample_variables(physeq)){
      message("The argument to `GroupBy` was:\n", GroupBy,
              "\nBut it was not found among sample variables:\n",
              paste0(sample_variables(physeq), collapse = ", "), "\n",
              "Please check the list shown above and try again.")
    }
  }
  # Start with fast melt
  mdt = fast_melt(physeq)
  if(!is.null(GroupBy)){
    # Add the variable indicated in `GroupBy`, if provided.
    sdt = data.table(SampleID = sample_names(physeq),
                     var1 = get_variable(physeq, GroupBy))
    setnames(sdt, "var1", GroupBy)
    # Join
    setkey(sdt, SampleID)
    setkey(mdt, SampleID)
    mdt <- sdt[mdt]
  }
  # Summarize
  nsamples=dim(otu_table(physeq))[2]
  summarydt = mdt[, list(meanRA = sum(RelativeAbundance)/nsamples,
                         sdRA = sd(RelativeAbundance),
                         minRA = min(RelativeAbundance),
                         maxRA = max(RelativeAbundance)),
                  by = c(Rank, GroupBy)]
  return(summarydt)
}

plot_taxa_summary = function(physeq, Rank, GroupBy = NULL){
  # Get taxa summary table 
  dt1 = summarize_taxa(physeq, Rank = Rank, GroupBy = GroupBy)
  # Set factor appropriately for plotting
  RankCol = which(colnames(dt1) == Rank)
  setorder(dt1, -meanRA)
  dt1[, RankFac := factor(dt1[[Rank]], 
                          levels = rev(dt1[[Rank]]))]
  dt1[, ebarMax := max(c(0, min(meanRA + sdRA))), by = eval(Rank)]
  dt1[, ebarMin := max(c(0, min(meanRA - sdRA))), by = eval(Rank)]
  # Set zeroes to one-tenth the smallest value
  ebarMinFloor = dt1[(ebarMin > 0), min(ebarMin)]
  ebarMinFloor <- ebarMinFloor / 10
  dt1[(ebarMin == 0), ebarMin := ebarMinFloor]
  
  pRank = ggplot(dt1, aes(x = meanRA, y = RankFac)) +
    scale_x_log10() +
    xlab("Mean Relative Abundance") +
    ylab(Rank) +
    theme_bw()
  if(!is.null(GroupBy)){
    # pRank <- pRank + facet_wrap(facets = as.formula(paste("~", GroupBy)))
    pRank <- pRank + geom_point(mapping = aes_string(colour = GroupBy),
                                size = 5)
  } else {
    # Don't include error bars for faceted version
    pRank <- pRank + geom_errorbarh(aes(xmax = ebarMax,
                                        xmin = ebarMin))
  }
  return(pRank)
}


plot_otu_on_map <- function(otu_column, map){
  colnames(otu_column)[4]="otu_abundance"
  otu_column$otu_abundance=as.numeric(as.character(otu_column$otu_abundance))
  colfunc<-colorRampPalette(c("dodgerblue2","khaki","orangered")) # create colours
  mycol <- function(x, myrange, n=100) round( 1+(x-myrange[1])/diff(myrange) * (n-1))
  ###whole range map#####
  #pdf("~/Dropbox/pathodopsis/pathodopsismicrobiome/figures/pathodopsis_collections_locations.pdf", width = 10 , height=8)
  map(database = "world",xlim = c(-15, 50), ylim = c(35, 65),col="#f0e4cd", fill= TRUE, resolution = 0, lty=1, lwd=0.2)
  map.axes()
  #("#003893","#fcd116","#ce1126","#0b8a46","#1ec5e4","#edbd39","#2162c8","#f3a782","#14b4c6","#145bc6","#14c67f","#c6145b","#767676","#d82129","#db977f","#1f174b","#dfa196","#ae030e","#8a6f7e","#d69031")
  ##############
  ###Cluster color is modified here, you can comment out next line if you don't want to customize your color i.e colors will be selected from the pallete above)
  #cluster_color_colors <- c("#CC0744", "#111111", "orangered", "#00846F", "#7B4F4B", "purple4","#A1C299" ,"deepskyblue4", '#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99','#b15928')
  cluster_color_colors <- c("#900c3f", "#204838", "#1f78b4", "#834848", "#505e5e", "#e31a1c", "#6a3d9a", "#b15928", "#ff7f00", "#35978f", "#d53e4f", "#bf812d", "#1a1a1a", "#4575b4", "#5aae61", "#dd3497", "#63531a", "#08306b", "#8c6bb1", "#d7301f")
  #points(x=test5$V3,y=test5$V2, pch=19,cex=1, col=(cluster_color_colors))
  #points(x=test5$V3,y=test5$V2, pch=as.integer(test5$V6),cex=1, col=(cluster_color_colors))
  
  #adjust transparency###
  colorList <- adjustcolor(cluster_color_colors[test5$V1], alpha.f=0.9)
  
  ###add points from all the sampled sites######

  ###Transparency#####
  colorGray <- adjustcolor("black", alpha.f=0.6)
  ###Points#####
  points(x=otu_column$Longitude  , y=otu_column$Latitude, cex=1, pch=20, bg=colorGray,  lwd=1, col=(otu_column$otu_abundance)) + scale_color_gradient()
  
  #dev.off()
  
  
}
