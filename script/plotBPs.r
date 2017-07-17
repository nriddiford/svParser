library(scales)
library(ggplot2)
library(dplyr)
library(RColorBrewer)


get_data <- function(infile = "all_bps_new.txt"){
  data<-read.delim(infile, header = F)
  colnames(data) <- c("event", "bp_no", "sample", "chrom", "bp", "gene", "feature", "type", "length")

  #filter on chroms
  data<-filter(data, chrom != "Y" & chrom != "4")
  #filter out samples
  data<-filter(data, sample != "A373R1" & sample != "A373R7" & sample != "A512R17" )
  data<-droplevels(data)
  dir.create(file.path("plots"), showWarnings = FALSE)
  return(data)
}


exclude_notch <- function(x){
  data<-get_data()
  data<-filter(data, !(chrom == "X" & bp >= 3000000 & bp <= 3300000))
  data<-droplevels(data)
  return(data)
}


clean_theme <- function(base_size = 12){
  theme(
    plot.title = element_text(hjust = 0.5, size = 20),
    panel.background = element_blank(),
    plot.background = element_rect(fill = "transparent",colour = NA),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.line.x = element_line(color="black", size = 0.5),
    axis.line.y = element_line(color="black", size = 0.5),
    axis.text = element_text(size=20),
    axis.title = element_text(size=30)
    )
}


set_cols <- function(df, col){
  names<-levels(df[[col]])
  cat("Setting colour levles:", names, "\n")
  level_number<-length(names)
  mycols<-brewer.pal(level_number, "Set2")
  names(mycols) <- names
  colScale <- scale_fill_manual(name = col,values = mycols)
  return(colScale)
}


plot_all_chroms_grid <- function(object=NA, notch=0){
  if(is.na(object)){
    object<-'type'
    cols<-set_cols(data, "type")
  }
  
  if(notch){
    data<-exclude_notch()
    ext<-'_excl.N.pdf'
  }
  else{
    data<-get_data()
    ext<-'.pdf'
  }

  cat("Plotting SVs by", object, "\n")
  
  p<-ggplot(data)
  p<-p + geom_histogram(aes(bp/1000000, fill = get(object)), binwidth=0.1, alpha = 0.8)  
  p<-p + facet_wrap(~chrom, scale = "free_x", ncol = 2)
  p<-p + scale_x_continuous("Mbs", breaks = seq(0,33,by=1), limits = c(0, 33),expand = c(0.01, 0.01))
  p<-p + scale_y_continuous("Number of Breakpoints", expand = c(0.01, 0.01))
  p<-p + clean_theme() +
    theme(axis.text.x = element_text(angle = 45, hjust=1),
          axis.text = element_text(size=12),
          axis.title = element_text(size=20),
          strip.text.x = element_text(size = 15)
    )
  
  if (object == 'type'){
    p<-p + cols
  }
  
  chrom_outfile<-paste("Breakpoints_chroms_by_", object, ext, sep = "")
  cat("Writing file", chrom_outfile, "\n")
  ggsave(paste("plots/", chrom_outfile, sep=""), width = 20, height = 10)
  
  p
}


bps_per_chrom <- function(object=NA, notch=0){
  if(is.na(object)){
    object<-'type'
    cols<-set_cols(data, "type")
  }
  
  if(notch){
    data<-exclude_notch()
    ext<-'_excl.N.pdf'
  }
  else{
    data<-get_data()
    ext<-'.pdf'
  }

  chromosomes<-c("2L", "2R", "3L", "3R", "X", "Y", "4")
  lengths<-c(23513712, 25286936, 28110227, 32079331, 23542271, 3667352, 1348131)

  karyotype<-setNames(as.list(lengths), chromosomes)

  for (c in chromosomes) {

    len<-karyotype[[c]]
    len<-len/1000000
    
    cat("Chrom", c, "length:", len, sep = " ",  "\n")

    per_chrom<-filter(data, chrom == c)

    p<-ggplot(per_chrom)
    p<-p + geom_histogram(aes(bp/1000000, fill = get(object)), binwidth=0.1, alpha = 0.8)
    p<-p + scale_x_continuous("Mbs", breaks = seq(0,len,by=1), limits = c(0, len+0.1),expand = c(0.01, 0.01))
    p<-p + scale_y_continuous("Number of Breakpoints", limits = c(0, 35), expand = c(0.01, 0.01))
    p<-p + clean_theme() +
      theme(axis.text.x = element_text(angle = 45, hjust=1),
        axis.title = element_text(size=20)
      )
    p<-p + ggtitle(paste("Chromosome: ", c))

    if (object == 'type'){
      p<-p + cols
    }
    
    per_chrom<-paste("Breakpoints_on_", c, "_by_", object, ext, sep = "")
    cat("Writing file", per_chrom, "\n")
    ggsave(paste("plots/", per_chrom, sep=""), width = 20, height = 10)
  }
}


bp_features <- function(notch=0){
  if(notch){
    data<-exclude_notch()
    ext<-'_excl.N.pdf'
  }
  else{
    data<-get_data()
    ext<-'.pdf'
  }
  
  # To condense exon counts into "exon"
  data$feature<-as.factor(gsub("_.*", "", data$feature))
  
  # Reoders descending
  data$feature<-factor(data$feature, levels = names(sort(table(data$feature), decreasing = TRUE)))
  
  cols<-set_cols(data, "feature")
  
  p<-ggplot(data)
  p<-p + geom_bar(aes(feature, fill = feature))
  p<-p + cols
  p<-p + clean_theme() +
    theme(axis.title.x=element_blank(),
      panel.grid.major.y = element_line(color="grey80", size = 0.5, linetype = "dotted"))
  p<-p + scale_x_discrete(expand = c(0.01, 0.01))
  p<-p + scale_y_continuous(expand = c(0.01, 0.01))

  features_outfile<-paste("Breakpoints_features_count", ext, sep = "")
  cat("Writing file", features_outfile, "\n")
  ggsave(paste("plots/", features_outfile, sep=""), width = 20, height = 10)

  p
}


sv_types<-function(notch=0,object=NA){
  if(is.na(object)){
    object<-'type'
  }
  
  if(notch){
    data<-exclude_notch()
    ext<-'_excl.N.pdf'
  }
  else{
    data<-get_data()
    ext<-'.pdf'
  }
  cols<-set_cols(data, "type")
  
  # Reorder by count
  data$type<-factor(data$type, levels = names(sort(table(data$type), decreasing = TRUE)))

  # Only take bp1 for each event
  data<-filter(data, bp_no != "bp2")

  p<-ggplot(data)
  p<-p + geom_bar(aes(get(object), fill = type))
  p<-p + cols
  p<-p + clean_theme() +
    theme(axis.title.x=element_blank(),
      panel.grid.major.y = element_line(color="grey80", size = 0.5, linetype = "dotted"),
	  axis.text.x = element_text(angle = 45, hjust=1),
	  axis.title = element_text(size=20)
    )
  p<-p + scale_x_discrete(expand = c(0.01, 0.01))
  p<-p + scale_y_continuous(expand = c(0.01, 0.01))

  types_outfile<-paste("sv_types_by_", object, ext, sep = "")
  cat("Writing file", types_outfile, "\n")
  ggsave(paste("plots/", types_outfile, sep=""), width = 20, height = 10)

  p
}


feature_lengths <- function(size_threshold = NA, notch=0){
  if(notch){
    data<-exclude_notch()
    ext<-'_excl.N.pdf'
  }
  else{
    data<-get_data()
    ext<-'.pdf'
  }
  
  cols<-set_cols(data, "type")
  
  # Only take bp1 for each event
  data<-filter(data, type != "TRA", type != "BND", bp_no != "bp2")

  data$length<-(data$length/1000)

  if(is.na(size_threshold)){
    size_threshold<-max(data$length)
  }
  
  ifelse(size_threshold <= 1, breaks<-0.1, breaks<-1)
  
  p<-ggplot(data, aes(length))
  p<-p + geom_density(aes(fill = type), alpha = 0.4)
  p<-p + clean_theme() +
    theme(panel.grid.major.y = element_line(color="grey80", size = 0.5, linetype = "dotted"))
  p<-p + scale_x_continuous("Size in Mb", expand = c(0,0), breaks = seq(0,size_threshold,by=breaks), limits=c(0, (size_threshold+0.1)))
  p<-p + scale_y_continuous(expand = c(0,0))
  p<-p + cols

  sv_classes_len_outfile <- paste("Classes_lengths", ext, sep = "")
  cat("Writing file", sv_classes_len_outfile, "\n")
  ggsave(paste("plots/", sv_classes_len_outfile, sep=""), width = 20, height = 10)

  p
}


feature_lengths_count <- function(size_threshold = NA, notch=0){
  if(notch){
    data<-exclude_notch()
    ext<-'_excl.N.pdf'
  }
  else{
    data<-get_data()
    ext<-'.pdf'
  }
  
  cols<-set_cols(data, "type")
  
  # Only take bp1 for each event
  data<-filter(data, type != "TRA", type != "BND", bp_no != "bp2")
  
  data$length<-(data$length/1000)
  
  if(is.na(size_threshold)){
    size_threshold<-max(data$length)
  }
  
  ifelse(size_threshold <= 1, breaks<-0.1, breaks<-1)
  
  p<-ggplot(data, aes(length))
  p<-p + geom_histogram(aes(length, ..count.., fill = type), colour = "black", binwidth = 0.05, position="dodge")
  p<-p + clean_theme() +
    theme(panel.grid.major.y = element_line(color="grey80", size = 0.5, linetype = "dotted"))
  p<-p + scale_x_continuous("Size in Mb", expand = c(0,0), breaks = seq(0,size_threshold,by=breaks), limits=c(0, (size_threshold+0.1)))
  p<-p + scale_y_continuous(expand = c(0,0))
  p<-p + geom_density(aes(fill = type),alpha=0.4, colour = NA)
  p<-p + cols
  
  sv_classes_len_outfile <- paste("Classes_lengths", ext, sep = "")
  cat("Writing file", sv_classes_len_outfile, "\n")
  ggsave(paste("plots/", sv_classes_len_outfile, sep=""), width = 20, height = 10)
  
  p
}


notch_hits <- function(){
  data<-get_data()
  data<-filter(data, chrom == "X", bp >= 3000000, bp <= 3300000)
  
  p<-ggplot(data)
  p<-p + geom_point(aes(bp/1000000, sample, colour = sample, shape = type, size = 2))
  p<-p + guides(color = FALSE, size = FALSE, sample = FALSE)
  p<-p + clean_theme() +
    theme(axis.title.y=element_blank(),
      panel.grid.major.y = element_line(color="grey80", size = 0.5, linetype = "dotted")
    )
  p<-p + scale_x_continuous("Mbs", expand = c(0,0), breaks = seq(3,3.3,by=0.05), limits=c(3, 3.301))

  p<-p + annotate("rect", xmin=3.000000, xmax=3.134532, ymin=0, ymax=0.5, alpha=.2, fill="green")
  p<-p + annotate("rect", xmin=3.134870, xmax=3.172221, ymin=0, ymax=0.5, alpha=.2, fill="skyblue")
  p<-p + annotate("rect", xmin=3.176440, xmax=3.300000, ymin=0, ymax=0.5, alpha=.2, fill="red")

  p
}


genome_hits <- function(notch=0){
  data<-get_data()
  
  if(notch){
    data<-exclude_notch()
    ext<-'_excl.N.pdf'
  }
  else{
    data<-get_data()
    ext<-'.pdf'
  }
  
  p<-ggplot(data)
  p<-p + geom_point(aes(bp/1000000, sample, colour = sample, shape = type, size = 0.5), alpha = 0.7)
  p<-p + guides(color = FALSE, size = FALSE)
  p<-p + clean_theme() +
    theme(axis.title.y=element_blank(),
      axis.text.x = element_text(angle = 45, hjust=1),
      axis.text = element_text(size=12),
      axis.title = element_text(size=20),
      strip.text.x = element_text(size = 15),
      panel.grid.major.y = element_line(color="grey80", size = 0.5, linetype = "dotted")
    )
  p<-p + facet_wrap(~chrom, scale = "free_x", ncol = 2)
  p<-p + scale_x_continuous("Mbs", breaks = seq(0,33,by=1), limits = c(0, 33),expand = c(0.01, 0.01))
  
  sv_gen_dist <- paste("bp_gen.dist", ext, sep = "")
  cat("Writing file", sv_gen_dist, "\n")
  ggsave(paste("plots/", sv_gen_dist, sep=""), width = 20, height = 10)
  
  p
}


feature_enrichment <- function(){
  genome_features<-read.delim('genome_features2.txt', header = T)
  data<-get_data()
  breakpoint_count<-nrow(data)
  
  # To condense exon counts into "exon"
  data$feature<-as.factor(gsub("_.*", "", data$feature))
  
  classes_count<-table(data$feature)
  
  class_lengths<-setNames(as.list(genome_features$length), genome_features$feature)
  cat("feature", "observed", "expected", "test", "sig", "p", "\n")
  
  for (f in levels(data$feature)) {
    feature_fraction<-class_lengths[[f]]/137547960
    feature_expect<-breakpoint_count*(class_lengths[[f]]/137547960)
    
    if(!is.null(class_lengths[[f]])){
      if(classes_count[f] >= feature_expect){
        stat<-binom.test(x = classes_count[f], n = breakpoint_count, p = feature_fraction, alternative = "greater")
        test<-"enrichment"
      }
      else{
        stat<-binom.test(x = classes_count[f], n = breakpoint_count, p = feature_fraction, alternative = "less")
        test<-"depletion"
      }
       
      ifelse(stat$p.value < 0.05, sig<-'T', sig<-'F')
       
      p_val<-format.pval(stat$p.value, digits = 3, eps=0.0001)
      cat(f, classes_count[f], feature_expect, test, sig, p_val, "\n")
      
      }
  } 
}

