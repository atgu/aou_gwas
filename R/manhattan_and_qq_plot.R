############## Load Packages ##############
packages <- c('ggplot2', 'tidyverse','dplyr','ggrepel','optparse','gridExtra','data.table','grDevices', 'stringr', 'purrr', 'readr', 'ggpubr')

for(p in packages){
  if(!require(p, character.only = T)){
    install.packages( p,  repos = c(CRAN = "http://cran.r-project.org") )
  }
}


############## Constants and util functions ##############
in2mm = 25.4
themes = theme_classic(base_size=8) + theme(plot.title = element_text(family = 'Arial', hjust = 0.5, color = 'Black', face = 'bold'),
               axis.text = element_text(family = 'Arial', color = 'Black'),
               axis.title = element_text(family = 'Arial', color = 'Black', face = 'bold'),
               legend.title = element_text(family = 'Arial', color = 'Black', face = 'bold'),
               legend.text = element_text(family = 'Arial', color = 'Black', size = 4),
               strip.text = element_text(family = 'Arial', color = 'Black', face = 'bold'),
               strip.background = element_rect( color = "black", size=0.5, linetype="solid") )

loglog_breaks = c(0:10, 15, 20, 50, 100, 200, 400)
ll_to_pvalue = function(yt) {
  return(if_else(yt >= 10, 10 ^ (yt / 10), yt))
}
pvalue_to_ll = function(y) {
  return(if_else(y >= 10, 10*log10(y), y))
}

gwas_loglog_trans <- function() {
  scales::trans_new("gwas_loglog", transform = pvalue_to_ll, inverse = ll_to_pvalue)
}


annotation_types = c('pLoF', 'missense|LC', 'synonymous', 'pLoF|missense|LC')
annotation_names = c('pLoF', 'Missense', 'Synonymous', 'pLoF|missense|LC')
names(annotation_names) = annotation_types
annotation_color_scale = scale_colour_manual(name = 'Annotation', values = colors, breaks = annotation_types, labels = annotation_names)
annotation_fill_scale = scale_fill_manual(name = 'Annotation', values = colors, breaks = annotation_types, labels = annotation_names)


format_chromosome <- function(chr_field){
  chr_field <- case_when(
    grepl('X', chr_field) ~ '23',
    grepl('Y', chr_field) ~ '24',
    grepl('M', chr_field) ~ '25',
    .default = str_replace(chr_field, 'chr', '')
  )
  chr_field <- as.numeric(chr_field)
  return(chr_field)
}

compute_lambda_with_af_cutoff <- function(data, pval_field, cutoff){
  lambda <- data %>%
    # dplyr::filter(get(af_field) >= cutoff) %>%
    dplyr::summarize(lambda = median(qchisq(get(pval_field), 1, lower.tail=F), na.rm=TRUE)/qchisq(0.5, 1))
  return(unlist(lambda$lambda))
}


############## Setup Arguments ##############
option_list <- list(
  make_option(c("-a","--alpha"), type="double", default=5e-8,
              help="defult significant level [default= %default]", metavar="double"),
  make_option(c("-b","--bp_col"), type="character",
              help="bp column [default= %default]", metavar="character"),
  make_option(c("-c","--chr_col"), type="character",
              help="chromosome column [default= %default]", metavar="character"),
  make_option(c("-d","--dn_cut"), type="double", default=0,
              help="Downsample cutoff [default= %default], -log10(p-value) below this cutoff value are down-sampled", metavar="double"),
  make_option(c("-e","--exp_pcol"), type="character", 
              default="Pvalue_expected_log10",
              # default=NULL,
              help="expected pvalue column [default= %default]", metavar="character"),
  make_option(c("-f", "--file"), type="character", default=NULL,
              help="dataset file name", metavar="character"),
  make_option(c("-g", "--gclambda"), type="double", default=NULL,
              help="LambdaGC", metavar="character"),
  make_option(c("-i","--locus"), type="character", default = 'locus',
              help="if given then chr:bp:ref:alt OR chr:bp identifier assumed and chr and bp are read from there [default= %default]", metavar="character"),
  make_option(c("-l", "--log10p"), type="logical", default=FALSE,
              help="whether the p-values are -log10 or not [default= %default]", metavar="logical"),
  make_option(c("-n", "--name"), type="character", default='gene_symbol',
              help="label field to print on manhattan plot [default= %default]", metavar="character"),
  make_option(c("-m", "--multi_pheno"), type="character",
              help="phenotype column [default= %default]", metavar="character"),
  make_option(c("-o", "--out"), type="character",
              help="output file name [default= %default]", metavar="character"),
  make_option(c("-p","--pval_col"), type="character", default="Pvalue_Burden_log10",
              help="pvalue column [default= %default]", metavar="character"),
  make_option(c("-q","--no_qqplot"), type="logical", default = FALSE,
              help="Making Manhattan plots only, ignore QQ plots [default= %default]", metavar="logical"),
  make_option(c("-s","--siglim"), type="integer", default=20,
              help="-upper limit for the number of significant bps to highlight [default= %default]", metavar="integer"),
  make_option(c("-t","--top"), type="integer", default=3,
              help="-number of top significant bps to highlight [default= %default]", metavar="integer"),
  make_option(c("-w","--af_col"), type="character", default='',
              help="allele count bin column [default= %default]", metavar="character")
);


############## Load Arguments ##############
opt_parser <- OptionParser(option_list=option_list);
opt <- parse_args(opt_parser, positional_arguments=0);

print(str(opt))
alpha <- opt$options$alpha
dn_cut <- opt$options$dn_cut
top <- opt$options$top
siglim <- opt$options$siglim
pcol <- opt$options$pval_col
no_qqplot <- opt$options$no_qqplot
af_col <- opt$options$af_col
log_p <- opt$options$log10p
locus_field <- opt$options$locus
name <- opt$options$name
exp_pcol <- opt$options$exp_pcol
lambda_gc <- opt$options$gclambda

file <- opt$options$file
print(paste("Reading file:", file))
data <- read_delim(file, delim = '\t')

if(!is.null(opt$options$bp_col) & !is.null(opt$options$chr_col)){
  bp_col <- opt$options$bp_col
  chr_col <- opt$options$chr_col
}else if(!is.null(opt$options$locus)){
  print("Getting chromosome and position from locus field")
  data <- data %>%
    dplyr::mutate(chr = str_split(get(locus_field), ':') %>% map_chr(., 1),
           pos = as.numeric(str_split(get(locus_field), ':') %>% map_chr(., 2)))
  bp_col <- 'pos'
  chr_col <- 'chr'
}else{
  stop('Locus information missing: provide input for (1) --locus or (2) --bp_col & --chr_col')
}

if(!is.null(opt$options$out)) {
  output_prefix <- opt$options$out
}else{
  output_prefix <- str_split(file, '\\.') %>% map_chr(., 1)
}

############## Chromosome Label Adjustment ##############
print('Formatting Chromosome column...')
data[[chr_col]] <- format_chromosome(data[[chr_col]])

############## Position Adjustment ##############
print('Formatting Base Pair Position column...')
data <- data[!is.na(data[[chr_col]]),]
data <- data[order(as.integer(data[[chr_col]])),]
chr <- unique(data[[chr_col]])
len <- vector()
ctr <- vector()
s <- 0

data[[bp_col]] <- as.numeric(data[[bp_col]])
print(head(data[[chr_col]]))
print(head(data[[bp_col]]))
for(i in chr){
  max <- max(data[data[[chr_col]] == i, bp_col])
  min <- min(data[data[[chr_col]] == i, bp_col])
  len[i] <- max
  data[data[[chr_col]] == i,"new_bp"] <- data[data[[chr_col]]== i, bp_col] + s
  s <- s + len[i]
  ctr[i] <- mean(c(max(data[data[[chr_col]] == i,'new_bp']),min(data[data[[chr_col]] == i,'new_bp'])))
}
names(ctr) <- gsub('24', 'Y', gsub('23', 'X',chr))

data$chr_col <- factor(data[[chr_col]],levels=chr)

if(af_col != ''){
  data$color_col <- data[[af_col]]
}else{
  data$color_col <- data$chr_col
}


############## Main function ##############
make_manhattan_and_qq_plot <- function(data, no_qqplot=no_qqplot, pcol=pcol, exp_pcol=exp_pcol, lambda_gc = lambda_gc,
                                       dn_cut=dn_cut, output_prefix=output_prefix, alpha=alpha, top=top, siglim=siglim, log_p=log_p, label=name){
  print(paste0('Formatting Pvalue column...'))
  data <- data %>%
      filter(!is.na(get(pcol)))
  n <- nrow(data)
  if(is.null(exp_pcol)){
     if(log_p){
       data[[pcol]] <- ifelse(10^-data[[pcol]] < 5e-324, 5e-324, 10^-data[[pcol]])
     }else{
       data[[pcol]] <- if_else(data[[pcol]]==0, 1e-323, data[[pcol]])
     }
    data$rank <- rank(data[[pcol]], ties.method = 'first')
    data <- data %>%
      dplyr::mutate(exp_pval = -log10(rank / (n+1)),
                    obs_pval = -log10(get(pcol)))
  }else{
    print(paste0('Obtain pre-computed expected Pvalues: ', exp_pcol, ', ', pcol))
    data$exp_pval <- data[[exp_pcol]]
    data$obs_pval <- data[[pcol]]
  }

  if(is.null(lambda_gc)){
    print(paste('Computing lambda GC'))
    lambda_gc <- median(qchisq(data[[pcol]], 1, lower.tail=F), na.rm=TRUE)/qchisq(0.5, 1)
  }


  run <- TRUE

  # Downsample observations with insignificant p-values
  if (dn_cut > 0){
    len <- sum(data$obs_pval <= dn_cut)
    print(paste0('Downsampling non-significant (p > 1e-', dn_cut ,') Pvalues (N = ', len, ')'))
    rm_ind <- order(data$obs_pval) %in% floor(seq.int(1,len,length.out = floor(max(0,len-1000))))
    subdata <- data[!rm_ind,] %>%
      mutate(label = factor(get(locus_field)))
    print(paste0('Raw total pvalues:', n))
    print(paste0('N pvalues removed:', length(rm_ind)))
    print(paste0('N pvalues kept:', nrow(subdata)))
  }else{
    subdata <- data
  }

  max <- if_else(is.infinite(max(subdata$obs_pval)), max(subdata[!is.infinite(subdata$obs_pval), 'obs_pval']), max(subdata$obs_pval))

  plt_qq <- subdata %>%
    # filter(!is.infinite(obs_pval)) %>%
    ggplot + aes(x=exp_pval, y=obs_pval) +
    geom_point(size = 0.1) +
    geom_abline(intercept = 0, slope = 1, lwd=0.15, color = '#CC3311')+
    labs(x=expression(bold(Expected -log[10](P))), y=expression(bold(Observed -log[10](P))))+
    themes +
    theme(panel.border = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank()) +
    annotate('text', x=0, y=max, label = (paste0("lambda[GC] : ", round(lambda_gc,digit=3))),
             parse = TRUE, color="#CC3311",vjust = 1, hjust = -0.1, size=2.5)
  if(!no_qqplot){
    png(paste0(output_prefix,"_",pcol,"_qqplot.png"), width=70/in2mm, height=70/in2mm, units = 'in', res = 300)
    print(paste0('Making QQ-plot:', output_prefix,"_",pcol,"_qqplot.png"))
    print(plt_qq)
    dev.off()
  }



  if('annotation' %in% colnames(subdata)){
    subdata <- subdata %>%
      mutate(annotation =factor(annotation, levels = annotation_types))
  }else{
    subdata <- subdata %>%
      filter(obs_pval>2)
  }
  plt_mh <- subdata %>%
    # filter(!is.infinite(obs_pval)) %>%
    ggplot + aes(x = new_bp, y = obs_pval, color = color_col, label= get(label)) +
    geom_point(size=0.1) +
    geom_hline(yintercept = -log10(alpha), lwd=0.15, color='#CC3311') +
    geom_hline(yintercept = 5, lwd=0.15, color='#e18a50') +
    scale_x_continuous(label = gsub('20|22', '',names(ctr)), breaks = ctr) +
    scale_y_continuous(trans = gwas_loglog_trans(), breaks=loglog_breaks)+
    scale_size_continuous(range = c(0.1,2)) +
    labs(x = "Chromosome",
         y = expression(bold(-log[10](P))))
  if(af_col == ''){
    plt_mh <- plt_mh +
      scale_color_manual(values =rep(c("#333366", "#6699CC"), length(chr))[1:length(chr)]) +
      theme(legend.position = 'None')
  }else{
    plt_mh <- plt_mh +
      scale_color_manual(name=NULL, values = c('#333366', '#6699CC', '#BFCFFF','#FFF381', '#F5A623', '#E56C4F')) +
      theme(legend.position = 'top', legend.direction = 'horizontal') +
      facet_wrap(~color_col, ncol = 1)
  }
  if(sum(subdata$obs_pval>-log10(alpha))==0){
    print("No significant value is presented")
    run <- FALSE
  } else if (sum(subdata$obs_pval>-log10(alpha))>top & sum(subdata$obs_pval>-log10(alpha))<siglim){
    subdata$top <- 0
    subdata$top[order(subdata$obs_pval, decreasing = T)[1:top]] <- 1
    subdata$sig <- as.integer(subdata$obs_pval>-log10(alpha))
  } else if (sum(subdata$obs_pval>-log10(alpha)) >= siglim) {
    subdata$top <- 0
    subdata$top[order(subdata$obs_pval, decreasing = T)[1:top]] <- 1
    subdata$sig <- 0
    subdata$sig[order(subdata$obs_pval, decreasing = T)[1:siglim]] <- 1
  } else if (sum(subdata$obs_pval>-log10(alpha)) > 0 & sum(subdata$obs_pval >-log10(alpha)) <= top ) {
    subdata$sig <- as.integer(subdata$obs_pval>-log10(alpha))
    subdata$top <- as.integer(subdata$obs_pval>-log10(alpha))
  }

  plt_mh <- plt_mh +
    themes

  if (run){
    subdata$sig2 <- subdata$sig==1 & subdata$top==0
    sub1 <- subset(subdata,sig2==1)
    sub2 <- subset(subdata,top==1)

    if (sum(subdata$sig)>top){
      plt_mh <- plt_mh +
        geom_point(data = sub1,color="#CC3311",alpha=0.7)+
        geom_point(data=sub2,color="#CC3311") +
        geom_label_repel(data = sub1,aes(fontface="bold"),color="#CC3311", size = 1.5 ,
                   min.segment.length = unit(0, 'lines'),
                   nudge_y = .2)+
        geom_label_repel(data = sub2,aes(fontface="bold"),color="#CC3311" , size = 2 ,
                   min.segment.length = unit(0, 'lines'),
                   nudge_y = .2)
    }else{
      plt_mh <- plt_mh +
        geom_point(data=sub2,color="#CC3311")+
        geom_label_repel(data = sub2,aes(fontface="bold"),color="#CC3311", size = 2)
    }
  }

  if(no_qqplot){
    height = if_else(af_col == '', 50, 250)
    png(paste0(output_prefix,"_",pcol,"_manhattan.png"), width=140/in2mm, height=height/in2mm, units = 'in', res = 300)
    print(plt_mh)
    dev.off()
  }else{
    figure = ggpubr::ggarrange(plt_mh, plt_qq, ncol = 2, widths = c(2.5, 1))
    png(paste0(output_prefix,"_",pcol,"_qq_and_manhattan.png"), width=174/in2mm, height=50/in2mm, units = 'in', res = 300)
    print(paste0('Combining the two figures: ', output_prefix,"_",pcol,"_qq_and_manhattan.png"))
    print(figure)
    dev.off()
  }
}




if(!is.null(opt$options$multi_pheno)){
  pheno_col <- opt$options$multi_pheno
  phenotypes <- names(table(data[[pheno_col]]))
  for(pheno in phenotypes){
    print(paste0('Phenotype: ', pheno))
    sub_data <- data %>%
      filter(get(pheno_col) == pheno)
    make_manhattan_and_qq_plot(data=sub_data, no_qqplot=no_qqplot, pcol=pcol,  dn_cut=dn_cut, output_prefix=paste0(output_prefix, '_', pheno), alpha=alpha, top=top, siglim=siglim, log_p=log_p, label=name)
  }
}else{
  if('annotation' %in% colnames(data)){
    for(annt in names(table(data$annotation))){
      subdata <- data%>% filter(annotation == annt)
      make_manhattan_and_qq_plot(data=subdata ,
                                 no_qqplot=no_qqplot,
                                 pcol=pcol,
                                 exp_pcol=exp_pcol,
                                 dn_cut=dn_cut,
                                 lambda_gc = lambda_gc,
                                 output_prefix=paste0(output_prefix, '_', annt),
                                 alpha=alpha,
                                 top=top,
                                 siglim=siglim,
                                 log_p=log_p,
                                 label=name
      )
    }
  }
  make_manhattan_and_qq_plot(data=data ,
                                 no_qqplot=no_qqplot,
                                 pcol=pcol,
                                 exp_pcol=exp_pcol,
                                 dn_cut=dn_cut,
                                 lambda_gc = lambda_gc,
                                 output_prefix=output_prefix,
                                 alpha=alpha,
                                 top=top,
                                 siglim=siglim,
                                 log_p=log_p,
                                 label=name)
}
