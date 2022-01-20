rm(list=ls())
library(openxlsx)
library(circlize)
library(ComplexHeatmap)
library(RColorBrewer)
library(ggplot2)

args = commandArgs(trailingOnly=TRUE)
if (length(args) < 2) {
    stop("Must provide at least two arguments: [1] Output directory, [2-...] TSVs from annotSV")
}
out_dir=args[1]
sample_files=args[2:length(args)]

#### FOR DRAGEN
min_span_reads=2   ## Minimum number of spanning reads for alt allele
min_split_reads=2   ## Mimimum number of split reads for alt allele


if (! dir.exists(out_dir)) {
  dir.create(out_dir, recursive = T)
}

#### Function to read AnnotSV results
read_annotsv <- function(annotfile,annotsv_suffix){
  
  #annotfile="/data/WES_breastcancer/analysis/annotSV/manta_results_annotated/Sample_5838+Sample_4333/Sample_5838+Sample_4333.somaticSV.annot.tsv"
  curr_sample=gsub(annotsv_suffix,"",basename(annotfile))
  # browser()
  print(paste0("opening ", annotfile))
  mydata <- read.table(file = annotfile, sep="\t", header=T, fill=F, quote="\"", stringsAsFactors = F)
  mydata$SV.size <- mydata$SV.end - mydata$SV.start
  
  sample_names_idx <- unique(c(which(colnames(mydata)=="FORMAT")+1, which(colnames(mydata)=="AnnotSV.type")-1))
  normal_tumor_ids <- colnames(mydata)[sample_names_idx]
  if (length(normal_tumor_ids) == 1) {
    names(normal_tumor_ids) <- c("tumor")
  }else {
    names(normal_tumor_ids) <- c("normal","tumor")
  }
  for (i in 1:length(normal_tumor_ids)) {
    id=normal_tumor_ids[i]
    split_fields <- lapply(strsplit(mydata[,id],":"), function(vec) {
      #return_vec = vec[5:6]
      return_vec = vec
      if (length(vec) == 1) {
        #print("Fixing one entry with no SR")
        return_vec <- c(vec, "0,0")
      }
      return(return_vec)
    })

    #pr_sr_counts <- data.frame(cbind(
    #  apply(do.call(rbind,strsplit(do.call(rbind,lapply(split_fields, "[[", 1)),",")),2,as.numeric),
    #  apply(do.call(rbind,strsplit(do.call(rbind,lapply(split_fields, "[[", 2)),",")),2,as.numeric)))
    pr_sr_counts <- data.frame(matrix(
          apply(do.call(rbind,strsplit(do.call(rbind,lapply(split_fields, "[[", 1)),",")),2,as.numeric),
          apply(do.call(rbind,strsplit(do.call(rbind,lapply(split_fields, "[[", 2)),",")),2,as.numeric),
      nrow=length(split_fields),ncol=4), stringsAsFactors=F)
    colnames(pr_sr_counts) <- paste0(names(normal_tumor_ids)[i], c("_pr_ref","_pr_alt","_sr_ref","_sr_alt"))
    # colnames(pr_sr_counts) <- c("spanningReads_ref","spanningReads_alt","splitReads_ref","splitReads_alt")

    mydata <- cbind(sample_id=id, mydata, pr_sr_counts)
    # mydata <- cbind(sample_id=id, mydata)
    colnames(mydata) <- gsub(id,"INFO",colnames(mydata))
  }

  return(mydata)
}


#### Read TSVs
all_tsv <- lapply(sample_files, read_annotsv, annotsv_suffix) 
all_manta_results.raw <- do.call(rbind, all_tsv)
all_manta_results <- all_manta_results.raw[all_manta_results.raw$AnnotSV.type=="full", ]

#### Filter AnnotSV results
filter_criteria <- all_manta_results$tumor_pr_alt >= min_span_reads | all_manta_results$tumor_sr_alt >= min_split_reads
all_manta_results <- all_manta_results[filter_criteria, ]
all_manta_results <- all_manta_results[all_manta_results$FILTER=="PASS", ]

if (nrow(all_manta_results) < 1) { stop("No SVs left after filtering.") }

library(openxlsx)
write.xlsx(all_manta_results, file=file.path(out_dir,"structural_variants.xlsx"))


#### Set up stuff to add gene name labels to the outside (experimental)
gene_info_cols <-c("SV.chrom","SV.start","SV.end","Gene.name")
genes_to_label <- unique(all_manta_results[,"Gene.name"])
#### Collapse gene labels if more than 3 genes
gene_info <- all_manta_results.raw[match(genes_to_label, all_manta_results.raw$Gene.name),gene_info_cols]
gene_info$gene_label <- unlist(lapply(gene_info$Gene.name, function(genelist) {
  split_gene=unlist(strsplit(genelist,"/"))
  if (length(split_gene) > 3) {
    label_str = paste0(length(split_gene), " genes")
  } else {
    label_str = paste0(split_gene, collapse="\n")
  }
  return(label_str)
}))

#### Define colors for the different SV types
all_sv_types <- sort(unique(all_manta_results$SV.type))
num_sv_types <- min(8,max(3, length(all_sv_types)))
sv_type_colors <- colorRampPalette(brewer.pal(num_sv_types, "Set1"))(num_sv_types)
names(sv_type_colors) <- sort(unique(all_manta_results$SV.type))

#### Make a simple barplot counting up the SVs by type
#svtypes <- as.data.frame.matrix(table(all_manta_results[,c("sample_id","SVTYPE")]))
svtypes <- as.data.frame.matrix(table(all_manta_results[,c("sample_id","SV.type")]))
svtypes$sample <- rownames(svtypes)
svtypes.melt <- reshape2::melt(svtypes, id.vars="sample", value.name="SV_Type")
svtypes.melt$sample <- factor(svtypes.melt$sample, levels=names(sort(table(all_manta_results$sample_id), decreasing = T)))
ggplot(svtypes.melt, aes(x=sample, y=SV_Type, fill=variable)) +
  xlab("")+ylab("Number of SVs") +
  scale_fill_manual(values=sv_type_colors) +
  geom_col() +
  theme_linedraw(base_size = 10) +
  theme(axis.text.x = element_text(angle=30, hjust=1))
ggsave(file.path(out_dir,"SV_counts.pdf"),height=5, width= min(c(50,length(unique(svtypes.melt$sample))*0.25)))


##########################################################################################
# This makes circos plots, but they don't look good/work for more than, say 6-8 samples...
##########################################################################################

###### Set up a bed-like data frame for input into circos.genomicTrackPlotRegion()
#bed <- data.frame(chr=paste0("chr",all_manta_results$SV.chrom),
#                 start=all_manta_results$SV.start,
#                 end=all_manta_results$SV.end,
#                 height=1,
#                 svtype=all_manta_results$SV.type,
#                 gene=all_manta_results$Gene.name,
#                 sample=all_manta_results$sample_id)
#
##### Split it into a list by sample (defines each concentric track within the circle)
#bed.split <- split(bed, bed$sample)
#
##### Make circos plot
##### Circle sections are chromosomes, each cell line in separate tracks within circle
##
#### First count up nubmer of SVs per sample per chromosome
#### Gonna use this to pick which chromosomes to show (showing top N most mutated chroms)
#### This is probably a stupid approach, but works for now
##top_n_chrs <- 5
##chr_count_tally <- as.data.frame.matrix(table(bed$sample, bed$chr))
## my_chrs <- names(tail(sort(colMeans(chr_count_tally),n = top_n_chrs)))
## my_chrs <- my_chrs[order(as.numeric(gsub("chr","",my_chrs)))]
##
#### Make a legend for the colors
##library(ComplexHeatmap)
#color_legend <- Legend(labels = names(sv_type_colors), type = "points",
#                      labels_gp = gpar(fontsize = 12),
#                      legend_gp = gpar(fill = sv_type_colors, col = sv_type_colors,pch=22),
#                      background = "white",
#                      # title_position = "topleft",
#                      title = "SV Type")
#
#### Make a bed-like object for gene labels from filtered chromosomes
##genecounts <- sort(table(unlist(lapply(bed.split, function(x){as.character(x$gene)}))))
##genes_to_label = names(genecounts[genecounts > (0.3*length(bed.split))])
##gene_bed <- gene_info[,c("SV.chrom","SV.start","SV.end","gene_label","Gene.name")]
##colnames(gene_bed)[1:3] <- c("chr","start","end")
##gene_bed$chr <- paste0("chr",gene_bed$chr)
##gene_bed <- gene_bed[gene_bed$Gene.name %in% genes_to_label &
##                       gene_bed$chr %in% my_chrs,]
##
#
#
#
##
################################################################################################
####################      PLOT #1: Top chromosomes + Gene Labels (testing)    ##################
################################################################################################
##pdf(paste0(out_dir,"/manta_circos.top.genelabels.pdf"),height=8, width=8)
##track_label_loc=-3e9  ## Location of the track labels; I haven't quite figured this out yet, but it seems to be a cumulative genomic coordinate (i.e. it sums up the genomic positions fromt the first chromosome it strats drawing)
##
#### Configure settings for margins, gaps, and canvas size before initializing
##circos.par(track.margin=c(0.0005,0.0005), 
##           gap.after = 5,
##           canvas.xlim = c(-1.2, 1.2),canvas.ylim = c(-1.2, 1.2))
##
#### Initialize the outer track with chromosomal ideograms
##circos.initializeWithIdeogram(cytoband="cytoBand.txt",
##                              plotType = "ideogram",
##                              track.height = 0.01,ideogram.height =0.025,
##                              chromosome.index = my_chrs)
#### Add gene labels (this is kinda crappy right now)
##circos.genomicLabels(gene_bed, labels.column = 4, side = "outside",cex=0.3)
##
#### This function makes blocks defined in the input (the bed-like data frame)
##myplotfun <- function(region, value, ...) {
##  circos.genomicRect(region, value[[1]],
##                     col = sv_type_colors[value[[2]]],
##                     border = sv_type_colors[value[[2]]],
##                     ...)
##}
##track_height=1/(length(bed.split)+4)  ### Defines the width of each concentric track
##                                      ### Need the extra padding (the "+4") to account for extra tracks like the ideogram and gene labels
##
#### Now, loop through the split bed list (i.e. each sample)
####  and add a rectangle/block for each SV
##for (i in 1:(length(bed.split))) {
##  curr_track=names(bed.split)[i]
##  require(circlize)
##  circos.genomicTrackPlotRegion(bed.split[[curr_track]], stack=T,
##                                track.height = track_height,
##                                bg.border="grey90",
##                                panel.fun = myplotfun)
##  ## This prints the name to the track
##  circos.text(track_label_loc,1,gsub("Sample.*_","",curr_track),
##              cex = 0.6,col="red3")
##}
##
###### Finally this adds the chromosome labels to the very outside
##circos.track(track.index = 1, 
##             panel.fun = function(x, y) {
##               circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
##                           facing = "outside", niceFacing = TRUE, adj = c(0.5, 2),
##                           cex=1.2, font=2)
##             }, bg.border = NA)
##
###### Render the legend
##draw(color_legend,x = unit(4, "mm"), y = unit(4, "mm"), just = c("left", "bottom"))
##
##title("Manta SV Calls")  ## Add a title
##circos.clear()   ## Not super necessary, but helpful when running small sections of the code
##dev.off()
##
##
#
#
#
#
#
###############################################################################################
###################          PLOT #2: Top chromosomes (No gene labels)       ##################
###############################################################################################
## pdf(paste0(out_dir,"/manta_circos.top.pdf"),height=8, width=8)
## track_label_loc=-3e9  ## Location of the track labels; I haven't quite figured this out yet, but it seems to be a cumulative genomic coordinate (i.e. it sums up the genomic positions fromt the first chromosome it strats drawing)
## 
## circos.par(track.margin=c(0.0005,0.0005), 
##            gap.after = 5,
##            canvas.xlim = c(-1.2, 1.2),canvas.ylim = c(-1.2, 1.2))
## circos.initializeWithIdeogram(cytoband="cytoBand.txt",
##                               plotType = "ideogram",
##                               track.height = 0.01,ideogram.height =0.025,
##                               chromosome.index = my_chrs)
## 
## myplotfun <- function(region, value, ...) {
##   circos.genomicRect(region, value[[1]],
##                      col = sv_type_colors[value[[2]]],
##                      border = sv_type_colors[value[[2]]],
##                      ...)
## }
## track_height=1/(length(bed.split)+2)
## 
## for (i in 1:(length(bed.split))) {
##   curr_track=names(bed.split)[i]
##   require(circlize)
##   circos.genomicTrackPlotRegion(bed.split[[curr_track]], stack=T,
##                                 track.height = track_height,
##                                 bg.border="grey90",
##                                 panel.fun = myplotfun)
##   circos.text(track_label_loc,1,gsub("Sample.*_","",curr_track),
##               cex = 0.6,col="red3")
## }
## circos.track(track.index = 1, 
##              panel.fun = function(x, y) {
##                circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
##                            facing = "outside", niceFacing = TRUE, adj = c(0.5, 2),
##                            cex=1.2, font=2)
##              }, bg.border = NA)
## 
## draw(color_legend,x = unit(4, "mm"), y = unit(4, "mm"), just = c("left", "bottom"))
## title("Manta SV Calls")
## circos.clear()
## dev.off()
## 
#
#
#
#
#
#
#
###############################################################################################
###################           PLOT #3: All chromosomes + No gene labels      ##################
###############################################################################################
#pdf(paste0(out_dir,"/manta_circos_all.pdf"),height=50, width=50)
#my_chrs <- sort(unique(bed$chr))
#
#track_label_loc=-1e9  ## Location of the track labels; I haven't quite figured this out yet, but it seems to be a cumulative genomic coordinate (i.e. it sums up the genomic positions fromt the first chromosome it strats drawing)
#circos.par(track.margin=c(0.0005,0.0005), 
#           gap.after = 5,
#           canvas.xlim = c(-1.2, 1.2),canvas.ylim = c(-1.2, 1.2))
#circos.initializeWithIdeogram(cytoband="cytoBand.txt",
#                              plotType = c("ideogram"),
#                              track.height = 0.01,ideogram.height =0.025,
#                              chromosome.index = my_chrs)
#
#myplotfun <- function(region, value, ...) {
#  circos.genomicRect(region, value[[1]],
#                     col = sv_type_colors[value[[2]]],
#                     border = sv_type_colors[value[[2]]],
#                     ...)
#  
#}
#track_height=1/(length(bed.split)+2)
#
#for (i in 1:length(bed.split)) {
#  curr_track=names(bed.split)[i]
#  require(circlize)
#  circos.genomicTrackPlotRegion(bed.split[[curr_track]], stack=T,
#                                track.height = track_height,
#                                bg.border="grey90",
#                                panel.fun = myplotfun)
#  # circos.text(track_label_loc,1,gsub("Sample.*_","",curr_track),
#  circos.text(track_label_loc,1,curr_track,
#              cex = 0.6,col="red3")
#}
#circos.track(track.index = 1, 
#             panel.fun = function(x, y) {
#               circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
#                           facing = "clockwise", niceFacing = TRUE, adj = c(-0.5, 0),
#                           cex=1.2, font=2)
#             }, bg.border = NA)
#draw(color_legend,x = unit(4, "mm"), y = unit(4, "mm"), just = c("left", "bottom"))
#title("Dragen SV Calls")
#circos.clear()
#dev.off()
#
#
#
