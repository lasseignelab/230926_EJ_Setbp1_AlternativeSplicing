### These functions are required for using the MARVEL vignette's way of
# using modified wiggleplotr functions to annotate splice junctions on full
# transcripts. MARVEL is available from CRAN, and wiggleplotr is available on
# bioconductor. These edited functions are not from the wiggleplotr source, they
# are from MARVEL's Sean Wen available at https://drive.google.com/file/d/1eEkq1axQXynVifWanh9ACichIrJjyNfr/view

## Emma Jones EDIT - enable non-wiggleplotr function names
adhocGene_PlotSJPosition_EJ <- function (MarvelObject, coord.intron, coord.intron.ext = 50, 
                                         rescale_introns = FALSE, show.protein.coding.only = TRUE,
                                         anno.label.size = 3, anno.colors = c("black", "gray", "red")) 
{
  MarvelObject <- MarvelObject
  coord.intron <- coord.intron
  sj.metadata <- MarvelObject$sj.metadata
  gtf <- MarvelObject$gtf
  coord.intron.ext <- coord.intron.ext
  rescale_introns <- rescale_introns
  show.protein.coding.only <- show.protein.coding.only
  anno.label.size <- anno.label.size
  anno.colors <- anno.colors
  gene_short_name <- sj.metadata[which(sj.metadata$coord.intron == 
                                         coord.intron), "gene_short_name.start"]
  if (length(gene_short_name) == 0) {
    message("No corresponding gene found for coord.intron specified")
    return(MarvelObject)
  }
  gtf <- gtf[grep(gene_short_name, gtf$V9, fixed = TRUE), ]
  message("Retrieving transcripts from GTF file...")
  . <- strsplit(gtf$V9, split = ";")
  . <- sapply(., function(x) grep("gene_name", x, value = TRUE))
  . <- gsub("gene_name", "", .)
  . <- gsub(" ", "", .)
  . <- gsub("\"", "", .)
  gtf$gene_short_name <- .
  gtf <- gtf[which(gtf$gene_short_name == gene_short_name), 
  ]
  . <- strsplit(gtf$V9, split = ";")
  . <- sapply(., function(x) grep("transcript_id", x, value = TRUE))
  . <- gsub("transcript_id", "", .)
  . <- gsub(" ", "", .)
  . <- gsub("\"", "", .)
  gtf$transcript_id <- .
  . <- strsplit(gtf$V9, split = ";")
  . <- sapply(., function(x) grep("transcript_biotype", x, 
                                  value = TRUE))
  . <- gsub("transcript_biotype", "", .)
  . <- gsub(" ", "", .)
  . <- gsub("\"", "", .)
  if (length(unique(.)) == 1 & unique(.)[1] == "character(0)") {
    . <- strsplit(gtf$V9, split = ";")
    . <- sapply(., function(x) grep("transcript_type", x, 
                                    value = TRUE))
    . <- gsub("transcript_type", "", .)
    . <- gsub(" ", "", .)
    . <- gsub("\"", "", .)
    gtf$transcript_biotype <- .
  }
  else {
    gtf$transcript_biotype <- .
  }
  . <- strsplit(gtf$V9, split = ";")
  . <- sapply(., function(x) grep("exon_id", x, value = TRUE))
  . <- gsub("exon_id", "", .)
  . <- gsub(" ", "", .)
  . <- gsub("\"", "", .)
  gtf$exon_id <- .
  anno <- unique(gtf[, c("transcript_id", "transcript_biotype")])
  anno <- anno[grep("ENST|ENSMUST", anno$transcript_id), ]
  message(paste(nrow(anno), " transcripts identified", sep = ""))
  transcript_ids <- anno$transcript_id
  grange.exon.list <- list()
  for (i in 1:length(transcript_ids)) {
    gtf.small <- gtf[which(gtf$transcript_id == transcript_ids[i]), 
    ]
    gtf.small <- gtf.small[which(gtf.small$V3 == "exon"), 
    ]
    grange <- GenomicRanges::GRanges(seqnames = S4Vectors::Rle(gtf.small$V1), 
                                     ranges = IRanges::IRanges(gtf.small$V4, width = (gtf.small$V5 - 
                                                                                        gtf.small$V4) + 1), strand = gtf.small$V7[1], 
                                     exon_id = gtf.small$exon_id, exon_name = gtf.small$exon_id, 
                                     exon_rank = c(1:length(gtf.small$exon_id)))
    grange.exon.list[[i]] <- grange
  }
  grange.exon.list <- GenomicRanges::GRangesList(grange.exon.list)
  names(grange.exon.list) <- transcript_ids
  grange.exon.sj.list <- list()
  transcript.ids <- NULL
  for (i in 1:length(grange.exon.list)) {
    grange <- grange.exon.list[[i]]
    exon <- as.data.frame(grange)
    exon$start <- as.numeric(exon$start)
    exon$end <- as.numeric(exon$end)
    . <- strsplit(coord.intron, split = ":", fixed = TRUE)[[1]]
    chr.sj <- .[1]
    start.sj <- as.numeric(.[2])
    end.sj <- as.numeric(.[3])
    exon.sj.start <- exon$end[which(exon$end == start.sj - 
                                      1)]
    exon.sj.end <- exon$start[which(exon$start == end.sj + 
                                      1)]
    exon.small <- exon[which(exon$end == exon.sj.start | 
                               exon$start == exon.sj.end), ]
    exon.small$start[which(exon.small$end == exon.sj.start)] <- exon.small$end[which(exon.small$end == 
                                                                                       exon.sj.start)] - coord.intron.ext
    exon.small$end[which(exon.small$start == exon.sj.end)] <- exon.small$start[which(exon.small$start == 
                                                                                       exon.sj.end)] + coord.intron.ext
    if (nrow(exon.small) != 0) {
      grange.sj <- GenomicRanges::GRanges(seqnames = S4Vectors::Rle(exon.small$seqnames), 
                                          ranges = IRanges::IRanges(exon.small$start, width = (exon.small$end - 
                                                                                                 exon.small$start) + 1), strand = exon.small$strand, 
                                          exon_id = exon.small$exon_id, exon_name = exon.small$exon_id, 
                                          exon_rank = c(1:length(exon.small$exon_id)))
      grange.exon.sj.list[[i]] <- grange.sj
      transcript.ids[i] <- names(grange.exon.list)[i]
    }
    else {
      grange.exon.sj.list[[i]] <- FALSE
      transcript.ids[i] <- FALSE
    }
  }
  index.keep <- which(transcript.ids != FALSE)
  grange.exon.sj.list <- grange.exon.sj.list[index.keep]
  transcript.ids <- transcript.ids[index.keep]
  if (length(grange.exon.sj.list) != 0) {
    grange.exon.sj.list <- GenomicRanges::GRangesList(grange.exon.sj.list)
    names(grange.exon.sj.list) <- transcript.ids
  }
  transcript_ids <- anno$transcript_id
  grange.cds.list <- list()
  transcript.ids <- NULL
  for (i in 1:length(transcript_ids)) {
    gtf.small <- gtf[which(gtf$transcript_id == transcript_ids[i]), 
    ]
    gtf.small <- gtf.small[which(gtf.small$V3 == "CDS"), 
    ]
    if (nrow(gtf.small) != 0) {
      grange <- GenomicRanges::GRanges(seqnames = S4Vectors::Rle(gtf.small$V1), 
                                       ranges = IRanges::IRanges(gtf.small$V4, width = (gtf.small$V5 - 
                                                                                          gtf.small$V4) + 1), strand = gtf.small$V7[1], 
                                       exon_id = gtf.small$exon_id, exon_name = gtf.small$exon_id, 
                                       exon_rank = c(1:length(gtf.small$exon_id)))
      grange.cds.list[[i]] <- grange
      transcript.ids[i] <- transcript_ids[i]
    }
    else {
      grange.cds.list[[i]] <- FALSE
      transcript.ids[i] <- FALSE
    }
  }
  index.keep <- which(transcript.ids != FALSE)
  grange.cds.list <- grange.cds.list[index.keep]
  transcript.ids <- transcript.ids[index.keep]
  if (length(grange.cds.list) != 0) {
    grange.cds.list <- GenomicRanges::GRangesList(grange.cds.list)
    names(grange.cds.list) <- transcript.ids
  }
  metadata <- data.frame(transcript_id = names(grange.exon.list), 
                         stringsAsFactors = FALSE)
  metadata$gene_short_name <- gene_short_name
  strand <- gtf$V7[1]
  if (strand == "+") {
    metadata$strand <- 1
  }
  else {
    metadata$strand <- -1
  }
  metadata <- join(metadata, anno, by = "transcript_id", type = "left")
  metadata$transcript_id.biotype <- paste(metadata$transcript_id, 
                                          " (", metadata$transcript_biotype, ")", sep = "")
  . <- data.frame(transcript_id = names(grange.exon.list), 
                  stringsAsFactors = FALSE)
  . <- join(., metadata[, c("transcript_id", "transcript_id.biotype")], 
            by = "transcript_id", type = "left")
  names(grange.exon.list) <- .$transcript_id.biotype
  if (length(grange.exon.sj.list) != 0) {
    . <- data.frame(transcript_id = names(grange.exon.sj.list), 
                    stringsAsFactors = FALSE)
    . <- join(., metadata[, c("transcript_id", "transcript_id.biotype")], 
              by = "transcript_id", type = "left")
    names(grange.exon.sj.list) <- paste(.$transcript_id.biotype, 
                                        "_SJ", sep = "")
  }
  . <- data.frame(transcript_id = names(grange.cds.list), stringsAsFactors = FALSE)
  . <- join(., metadata[, c("transcript_id", "transcript_id.biotype")], 
            by = "transcript_id", type = "left")
  names(grange.cds.list) <- .$transcript_id.biotype
  metadata$transcript_id <- metadata$transcript_id.biotype
  metadata$transcript_id.biotype <- NULL
  metadata$transcript_biotype <- NULL
  if (length(grange.exon.sj.list) != 0) {
    grange.exon.list <- c(grange.exon.list, grange.exon.sj.list)
  }
  if (length(grange.exon.sj.list) != 0) {
    metadata <- data.frame(transcript_id = names(grange.exon.list), 
                           gene_short_name = gene_short_name, stringsAsFactors = FALSE)
    strand <- gtf$V7[1]
    if (strand == "+") {
      metadata$strand <- 1
    }
    else {
      metadata$strand <- -1
    }
  }
  if (show.protein.coding.only == TRUE) {
    transcript_ids <- metadata[grep("protein_coding", metadata$transcript_id, 
                                    fixed = TRUE), "transcript_id"]
    if (length(transcript_ids) != 0) {
      metadata <- metadata[grep("protein_coding", metadata$transcript_id, 
                                fixed = TRUE), ]
      grange.exon.list <- grange.exon.list[metadata$transcript_id]
      overlap <- intersect(names(grange.cds.list), transcript_ids)
      grange.cds.list <- grange.cds.list[overlap]
    }
    else {
      message("No protein-coding transcripts found for this gene")
      MarvelObject$adhocGene$SJPosition$metadata <- metadata
      return(MarvelObject)
    }
  }
  plot <- plotTranscripts(exons = grange.exon.list, 
                          cdss = grange.cds.list, transcript_annotations = metadata, 
                          rescale_introns = rescale_introns, new_intron_length = 50, 
                          flanking_length = c(50, 50), connect_exons = TRUE, transcript_label = TRUE, 
                          region_coords = NULL, anno.colors = anno.colors, anno.label.size = anno.label.size)
  MarvelObject$adhocGene$SJPosition$Plot <- plot
  MarvelObject$adhocGene$SJPosition$metadata <- metadata
  MarvelObject$adhocGene$SJPosition$exonfile <- grange.exon.list
  MarvelObject$adhocGene$SJPosition$cdsfile <- grange.cds.list
  return(MarvelObject)
}

## Use raw code for custom function

plotTranscripts <- function(exons, cdss = NULL, transcript_annotations = NULL, 
                            rescale_introns = TRUE, new_intron_length = 50, 
                            flanking_length = c(50,50), connect_exons = TRUE, 
                            transcript_label = TRUE, region_coords = NULL, anno.colors=NULL, anno.label.size=3){
  
  # Example arguments
  #library(wiggleplotr)
  #library(dplyr)
  #library(GenomicRanges)
  #library(GenomicFeatures)
  #library(biomaRt)
  #library(ggplot2)
  #library(ggpattern)
  
  #source("/Users/seanwen/Documents/MARVEL/wiggleplotr/wiggleplotr/R/functions.R")
  #source("/Users/seanwen/Documents/MARVEL/wiggleplotr/wiggleplotr/R/shortenIntrons.R")
  #source("/Users/seanwen/Documents/MARVEL/wiggleplotr/wiggleplotr/R/makePlots.R")
  
  #exons <- grange.exon.list
  #cdss <- grange.cds.list
  #transcript_annotations <- metadata
  #rescale_introns <- TRUE
  #new_intron_length <- 50
  #flanking_length <- c(50,50)
  #connect_exons <- TRUE
  #transcript_label <- TRUE
  #region_coords <- NULL
  #anno.colors <- c("black", "grey", "red")
  #anno.label.size <- 3
  
  ######################################################################
  
  #IF cdss is not specified then use exons instead on cdss
  if(is.null(cdss)){
    cdss = exons
  }
  
  #Check exons and cdss
  assertthat::assert_that(is.list(exons)|| is(exons, "GRangesList")) #Check that exons and cdss objects are lists
  assertthat::assert_that(is.list(cdss) || is(exons, "GRangesList"))
  
  #Join exons together
  joint_exons = joinExons(exons)
  
  #Extract chromosome name
  chromosome_name = as.vector(GenomicRanges::seqnames(joint_exons)[1])
  
  #If region_coords is specificed, then ignore the flanking_length attrbute and compute
  # flanking_length form region_coords
  if(!is.null(region_coords)){
    gene_range = constructGeneRange(joint_exons, c(0,0))
    min_start = min(GenomicRanges::start(gene_range))
    max_end = max(GenomicRanges::end(gene_range))
    flanking_length = c(min_start - region_coords[1], region_coords[2] - max_end)
  }
  #Make sure that flanking_length is a vector of two elements
  assertthat::assert_that(length(flanking_length) == 2) 
  
  #Rescale introns
  if (rescale_introns){
    tx_annotations = rescaleIntrons(exons, cdss, joint_exons, new_intron_length = new_intron_length, flanking_length)
    xlabel = "Distance from region start (bp)"
  } else {
    old_introns = intronsFromJointExonRanges(GenomicRanges::ranges(joint_exons), flanking_length = flanking_length)
    tx_annotations = list(exon_ranges = lapply(exons, GenomicRanges::ranges), cds_ranges = lapply(cdss, GenomicRanges::ranges),
                          old_introns = old_introns, new_introns = old_introns)
    
    xlabel = paste("Chromosome", chromosome_name, "position (bp)")
  }
  
  #If transcript annotations are not supplied then construct them manually from the GRanges list
  if(is.null(transcript_annotations)){
    plotting_annotations = dplyr::tibble(transcript_id = names(exons),
                                         strand = extractStrandsFromGrangesList(exons)) %>%
      prepareTranscriptAnnotations()
  } else{
    plotting_annotations = prepareTranscriptAnnotations(transcript_annotations)
  }
  
  #Plot transcript structures
  limits = c( min(IRanges::start(tx_annotations$new_introns)), max(IRanges::end(tx_annotations$new_introns)))
  structure = prepareTranscriptStructureForPlotting(tx_annotations$exon_ranges, 
                                                    tx_annotations$cds_ranges, plotting_annotations)
  
  ################################### SEAN EDIT ###################################
  
  # Indicate SJ
  structure$feature_type[grep("SJ", structure$transcript_label)] <- "sj"
  levels <- intersect(c("exon", "cds", "sj"), unique(structure$feature_type))
  structure$feature_type <- factor(structure$feature_type, levels=levels)
  
  # Create transcript id column without SJ suffix
  structure$transcript_id. <- gsub("_SJ", "", structure$transcript_id, fixed=TRUE)
  
  # Reorder transcripts, SJ: So that SJ is above transcript
  structure$transcript_id. <- factor(structure$transcript_id., levels=unique(structure$transcript_id.))
  
  transcript.id.sj <- unique(structure[which(structure$feature_type=="sj"), "transcript_id."])
  
  if(length(transcript.id.sj) != 0) {
    
    # Split data frame into transcripts with SJ, without SJ
    structure. <- structure[which(structure$transcript_id. %in% transcript.id.sj),]
    structure <- structure[-which(structure$transcript_id. %in% transcript.id.sj),]
    
    .list <- list()
    
    for(i in 1:length(transcript.id.sj)) {
      
      structure.small <- structure.[which(structure.$transcript_id.==transcript.id.sj[i]), ]
      
      # Re-label transcripts
      transcript.label <- unique(grep("_SJ", structure.small$transcript_label, fixed=TRUE, value=TRUE, invert=TRUE))
      structure.small$transcript_label[grep("_SJ", structure.small$transcript_label, fixed=TRUE, invert=TRUE)] <- ""
      structure.small$transcript_label[grep("_SJ", structure.small$transcript_label, fixed=TRUE)] <- transcript.label
      
      # Close gap (rank) between SJ and transcript
      structure.small$transcript_rank[which(structure.small$transcript_rank==max(structure.small$transcript_rank))] <- structure.small$transcript_rank[which(structure.small$transcript_rank==max(structure.small$transcript_rank))] - 0.5
      
      # Save into list
      .list[[i]] <- structure.small
      
    }
    
    # Merge
    structure. <- do.call(rbind.data.frame, .list)
    structure <- rbind.data.frame(structure., structure)
    
    # Standardise transcript labels positions
    #strand <- ifelse(structure$strand[1]=="-1", "-", "+")
    
    #if(strand=="-") {
    
    structure$label.start <- min(structure$start, na.rm=TRUE)
    
    #} else {
    
    #structure$label.start <- max(structure$start, na.rm=TRUE)
    
    #}
    
    
  }
  
  # Remove intermediate columns
  structure$transcript_id. <- NULL
  
  ################################### SEAN EDIT ###################################
  
  if(is.null(anno.colors[1])) {
    
    gg_color_hue <- function(n) {
      hues = seq(15, 375, length = n + 1)
      hcl(h = hues, l = 65, c = 100)[1:n]
    }
    
    n = length(levels(structure$feature_type))
    anno.colors = gg_color_hue(n)
    
  } else {
    
    anno.colors <- anno.colors
    
  }
  
  #source("/Users/seanwen/Documents/MARVEL/wiggleplotr/wiggleplotr/R/makePlots.R")
  plot = plotTranscriptStructure(structure, limits, connect_exons = connect_exons, xlabel = xlabel,
                                 transcript_label = transcript_label, anno.colors=anno.colors, anno.label.size=anno.label.size)
  return(plot)
}

## Get joinExons function
joinExons <- function(exons) {
  #Join a list of exons into one GRanges object
  
  #Test that all transcripts are on the same chromosome
  chrs = purrr::map_chr(as.list(exons), ~GenomicRanges::seqnames(.)[1] %>% 
                          S4Vectors::as.vector.Rle(mode = "character"))
  if (!all(chrs == chrs[1])){
    stop("Some transcripts are on different chromosomes.")
  }
  
  #Join all exons together
  transcript_ids = names(exons)
  joint_exons = c()
  for(tx_id in transcript_ids){
    tx = exons[[tx_id]]
    if(length(joint_exons) == 0){
      joint_exons = tx
    }
    else{
      joint_exons = c(joint_exons, tx)
    }
  }
  joint_exons = GenomicRanges::reduce(joint_exons)
  return(joint_exons)
}

## Get intronsFromJointExonRanges function

intronsFromJointExonRanges <- function(joint_exon_ranges, flanking_length){
  #Construct intron ranges from joint exon ranges
  introns = IRanges::gaps(joint_exon_ranges, 
                          start = min(IRanges::start(joint_exon_ranges)) - flanking_length[1], 
                          end = max(IRanges::end(joint_exon_ranges)) + flanking_length[2])
  return(introns)
}

## Get prepareTranscriptAnnotations

prepareTranscriptAnnotations <- function(transcript_annotations){
  assertthat::assert_that(assertthat::has_name(transcript_annotations, "transcript_id"))
  assertthat::assert_that(assertthat::has_name(transcript_annotations, "strand"))
  
  
  #Make sure that the strand information is represented correctly
  transcript_annotations = dplyr::mutate(transcript_annotations,
                                         strand = ifelse(strand %in% c("+","*") | strand == 1, 1, -1))
  
  #Add transcript label
  if(assertthat::has_name(transcript_annotations, "gene_name")){
    transcript_annotations = dplyr::select_(transcript_annotations, "transcript_id", "gene_name", "strand") %>% 
      dplyr::mutate(transcript_label = ifelse(strand == 1, 
                                              paste(paste(gene_name, transcript_id, sep = ":")," >",sep =""), 
                                              paste("< ",paste(gene_name, transcript_id, sep = ":"),sep ="")))
  } else{
    transcript_annotations = dplyr::mutate(transcript_annotations, transcript_label = ifelse(strand == 1, 
                                                                                             paste(paste(transcript_id, sep = ":")," >",sep =""), 
                                                                                             paste("< ",paste(transcript_id, sep = ":"),sep =""))) 
  }
  return(transcript_annotations)
}

## Get prepareTranscriptStructureForPlotting

prepareTranscriptStructureForPlotting <- function(exon_ranges, cds_ranges, transcript_annotations){
  #Combine exon_ranges and cds_ranges into a single data.frame that also contains transcript rank
  
  #Convert exon ranges into data.frame and add transcript rank
  exons_df = purrr::map_df(exon_ranges, data.frame, .id = "transcript_id")
  exons_df = dplyr::mutate(exons_df, transcript_rank = as.numeric(factor(exons_df$transcript_id)), type = "")
  transcript_rank = unique(exons_df[,c("transcript_id", "transcript_rank", "type")])
  
  #Convert CDS ranges into a data.frame
  cds_df = purrr::map_df(cds_ranges, data.frame, .id = "transcript_id")
  cds_df = dplyr::left_join(cds_df, transcript_rank, by = "transcript_id") #Add matching transcript rank
  
  #Join exons and cdss together
  exons_df = dplyr::mutate(exons_df, feature_type = "exon")
  cds_df = dplyr::mutate(cds_df, feature_type = "cds")
  transcript_struct = rbind(exons_df, cds_df)
  
  #Add transcript label to transcript structure
  transcript_struct = dplyr::left_join(transcript_struct, transcript_annotations, by = "transcript_id")
  return(transcript_struct)
}

## Get plotTranscriptStructure

plotTranscriptStructure <- function(exons_df, limits = NA, connect_exons = TRUE,  
                                    xlabel = "Distance from gene start (bp)", transcript_label = TRUE, anno.colors=NULL, anno.label.size=2){
  
  
  # Example arguments (Sean Edit)
  #exons_df <- structure
  #connect_exons <- TRUE
  #xlabel = "Distance from gene start (bp)"
  #transcript_label = TRUE
  #anno.colors <- anno.colors
  #anno.label.size <- 2
  
  #Extract the position for plotting transcript name
  transcript_annot = dplyr::group_by_(exons_df, ~transcript_id) %>% 
    dplyr::filter_(~feature_type %in% c("exon", "sj")) %>%
    dplyr::arrange_('transcript_id', 'start') %>%
    dplyr::filter(row_number() == 1)
  
  # Color scheme (Sean edit)
  if(is.null(anno.colors[1])) {
    
    gg_color_hue <- function(n) {
      hues = seq(15, 375, length = n + 1)
      hcl(h = hues, l = 65, c = 100)[1:n]
    }
    
    n = length(levels(structure$feature_type))
    anno.colors = gg_color_hue(n)
    
  } else {
    
    anno.colors <- anno.colors
    
  }
  
  #Create a plot of transcript structure
  plot = ggplot(exons_df) + geom_blank()
  if(connect_exons){ #Print line connecting exons
    plot = plot + geom_line(aes_(x = ~start, y = ~transcript_rank, group = ~transcript_rank, color = ~feature_type))
  }
  plot = plot + 
    geom_rect(aes_(xmin = ~start, 
                   xmax = ~end, 
                   ymax = ~transcript_rank + 0.25, 
                   ymin = ~transcript_rank - 0.25, 
                   fill = ~feature_type)) +
    theme_light() +
    theme(plot.margin=unit(c(0,1,1,1),"line"), 
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          legend.position="none",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.text.y = element_text(colour = "grey10"),
          strip.background = element_rect(fill = "grey85")) +
    xlab(xlabel) +
    facet_grid(type~.) +
    scale_y_continuous(expand = c(0.2,0.15)) +
    scale_fill_manual(values = anno.colors) +
    scale_colour_manual(values = anno.colors)
  if(all(!is.na(limits))){
    plot = plot + scale_x_continuous(expand = c(0,0)) +
      coord_cartesian(xlim = limits)
  }
  if(transcript_label){
    plot = plot + geom_text(aes_(x = ~label.start,
                                 y = ~transcript_rank + 0.30, 
                                 label = ~transcript_label),
                            data = transcript_annot, hjust = 0, vjust = 0, size = anno.label.size)
    
  }
  return(plot)
}
