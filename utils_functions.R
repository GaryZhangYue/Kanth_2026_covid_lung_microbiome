# Functions

# function to set a column as rownames and remove that column from the df; default is setting the first column as the rowname
setrowname =  function(df,col = 1){
  if (is.numeric(col)) {
    # print('it is a number')
    rownames(df) = df[,col]
    df = df[,-col]
  } else {
    # print('it is a string')
    rownames(df) = df[,which(names(df) %in% c(col))]
    df = df[,-which(names(df) %in% c(col))]
  }
  return(df)
} 


# function to generate collapsed taxa bar plot
plot_collapsed_taxabar = function(taxatable = asv.taxa.phylum[[1]],grp,rank = 'Phylum', ft = 9){
  # filtered metadata
  metatable.tmp = subset(metadata.o,get(grp) != 'N.A.', select = grp)
  # filtered taxa table
  asv.offp = select(taxatable,rownames(metatable.tmp))
  stopifnot(all(round(colSums(asv.offp)*100) == 100)) # sanity check
  
  #retain the ft most abundant taxa; calculate the sum of all the others
  asv.offp = asv.offp[order(rowSums(asv.offp), decreasing = T),]  
  asv.offps = bind_rows(asv.offp[c(1:ft),],colSums(asv.offp[c(c(ft+1):nrow(asv.offp)),]))
  rownames(asv.offps)[c(ft+1)] = 'All others'
  # stack it
  if (all(round(colSums(asv.offps)*100) ==100)) {
    rownames(asv.offps) -> asv.offps[[rank]]
    asv.offps.s = gather(asv.offps,"SD","frequency",-{{rank}})
    #merge the stacked ASV table with a subset of metadata to include some categorization columns
    metatable.tmp %>%
      merge(x = ., y = asv.offps.s, all.x = T, by.x = 0, by.y = "SD") -> asv.offps.s
    #coerce the phylum column to ordered factors for plotting
    asv.offps.s[[rank]] = factor(asv.offps.s[[rank]], levels = c(rownames(asv.offps)))
    asv.offps.s$Row.names = as.character(asv.offps.s$Row.names)
  } else {asv.offps.s = NULL}
  
  # collapse to groups
  asv.offps.s.collapsed = asv.offps.s %>%
    group_by_at(c(grp,rank) ) %>%
    summarise(frequency_avg=mean(frequency))
  
  p = ggplot(data = asv.offps.s.collapsed, mapping = aes(x = get(grp), y = `frequency_avg`,fill = get(rank))) +
    geom_bar(stat="identity",color = 'black') +
    scale_fill_manual(values = col21) +
    theme_classic() +
    theme(legend.title = element_blank(),
          legend.position = 'right',
          legend.text = element_text(size = 12),
          legend.key.size = unit(0.3,'cm'),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=10, face = 'bold', color = 'black'),
          axis.text.y = element_text(size=15, face = 'bold', color = 'black'),
          axis.title=element_text(size=18,face="bold"),
          strip.text.x = element_text(size = 15, color = "black", face = "bold"),
          strip.text.y = element_text(size = 15, color = "black", face = "bold"),
          strip.background = element_rect(color="black", fill="grey95", size=1.5, linetype="solid"))+
    xlab("Groups") + ylab('Relative Abundance') + labs(title = grp) +
    guides(fill=guide_legend(ncol=1))
  print(p)
  ggsave(plot = p,paste(PFIG,"/",'Taxonomy_top_9_phyla_barplot_', grp, '.pdf', sep = ""),width = 6, height = 4, device = 'pdf')
  
  return(list(asv.offps,asv.offps.s.collapsed))
}


# define a function for easy plotting
plot_alpha = function(altab = alpha, grp,index){
  altab = subset(altab, get(grp) != 'N.A.', select = c(grp,index))
  # Calculate dynamic ylim
  ylimitmax <- ceiling(range(altab[[index]])[2]*1.2)
  label_y_position <- ylimitmax*0.99  # Position label just below the max value
  
  plt = ggboxplot(altab, x = grp, y = index, fill = grp,  # change here 3, obj, y = ...
                  size = 1, add = "jitter")+
    xlab("Groups") + ylab(index) + ggtitle(paste0(index, ': ', grp)) + # change here ylab and ggtitle
    coord_cartesian(ylim = c(0,ylimitmax))+ # change here max()
    theme(text = element_text(size = 15), legend.position = 'none',
          axis.text.x = element_text(angle = 15, hjust = 1))+
    scale_fill_lancet(alpha=0.8)
  
  # run wilcox for comparing 2 groups and run kruskal wallis test for comparing more than 2 groups
  if (length(unique(altab[[grp]]))>2) {
    plt = plt + stat_compare_means(method = "kruskal.test",label.x = 1,label.y = label_y_position,size = 6)
  } else {
    plt = plt + stat_compare_means(method = "wilcox.test",label.x = 1,label.y = label_y_position,size = 6)
  }
  print(plt)
}


# function for easy plotting pcoa plots (continuous options might not work)
plot_unifrac = function(dm,metatable,grp,dm_name,permtest = F,grp_continuous = F,custom_colors){
  set.seed(666)
  # dummy input variable for testing function
  # dm = distmtx
  # metatable = metadata.o
  # dm_name = 'x'
  # grp = grouping
  
  
  metatable = subset(metatable, get(grp) != 'N.A.', select = c(grp))
  # select distance matrix matching to metadata
  dm = dist_subset(dm,rownames(metatable))
  # some distance measures may result in negative eigenvalues. In that case, add a correction:
  pcoa.wunifrac = pcoa(dm,correction = 'cailliez')
  # transform the pcoa.wunifrac object into a df for making it plottable with ggplot2
  coordinates.pcoa.wunifrac = data.frame(pcoa.wunifrac$vectors,check.names = FALSE)
  
  # Adonis test (only run if all groups have >= 5 samples)
  ## Calculate the number of samples in each group
  group_counts <- table(metatable[[grp]])
  ## Check if any group has fewer than 5 samples
  any_small_groups <- any(group_counts < 5)
  if (any_small_groups == F & permtest ==T){
    stopifnot(all(rownames(metatable) == names(dm)))
    adonis.r = adonis2(dm~get(grp),metatable,permut = 1000,parallel = 1)
    #adonis.p = round(adonis.r$`Pr(>F)`[1],digit = 4)
    adonis.p = adonis.r$`Pr(>F)`[1]
    
    betadisp.r = permutest(betadisper(dm,metatable[[grp]]),permutations = 1000,parallel = 1)
    #betadisp.p = round(betadisp.r$tab[[6]][1],digit = 4)
    betadisp.p = betadisp.r$tab[[6]][1]
    
    subtitle = paste0('permanova p val = ', adonis.p, ',betadisp p val = ', betadisp.p)
  } else {
    subtitle = 'Adonis/permdisp skipped'
  }
  
  
  # plot
  if (!all(rownames(metatable) == rownames(coordinates.pcoa.wunifrac))) {stop("mismatchs in pca coordinates and coloring column")} # sanity check
  
  # Calculate group centroids
  coordinates = coordinates.pcoa.wunifrac[,c(1,2)]
  #coordinates[[grp]]= metatable[[grp]]      
  coordinates[[grp]] = metatable[[grp]][match(rownames(coordinates),rownames(metatable))] 
  #print(coordinates)
  centroids <- coordinates %>%
    group_by_at(grp) %>%
    summarize(CentroidX = mean(Axis.1), CentroidY = mean(Axis.2))
  coordinates = merge(coordinates,centroids,
                      by = grp)           
  # automatically grab relative eigenvalue
  if ("Rel_corr_eig" %in% colnames(pcoa.wunifrac$values)) {
    pc_exp <- pcoa.wunifrac$values$Rel_corr_eig
  } else {
    pc_exp <- pcoa.wunifrac$values$Relative_eig
  }
  
  if (grp_continuous == F){
    
    pcoa <- ggplot() +
      # Smaller scatter points
      geom_point(data = coordinates,
                 mapping = aes(x = `Axis.1`, y = `Axis.2`, colour = get(grp)),
                 alpha = 0.6, size = 0.5) + 
      
      # Connection lines to centroids
      geom_segment(data = coordinates,
                   aes(x = CentroidX, xend = Axis.1, y = CentroidY, yend = Axis.2, color = get(grp)),
                   alpha = 0.3, size = 0.3) +
      # Add ellipses
      #stat_ellipse(data = coordinates,
      #              aes(x = `Axis.1`, y = `Axis.2`, color = get(grp)),
      #               type = "norm", level = 0.68, size = 0.5, alpha = 0.4) +
      # Emphasize centroids
      geom_point(data = centroids,
                 aes(x = CentroidX, y = CentroidY, fill = get(grp)),
                 size = 4, shape = 21, stroke = 0.6, color = 'black') +
      
      labs(title = paste0(dm_name, ': ', grp),
           subtitle = subtitle) +
      
      xlab(paste('PC1 (', label_percent()(pc_exp[1]), ')', sep = '')) +
      ylab(paste('PC2 (', label_percent()(pc_exp[2]), ')', sep = '')) +
      
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            axis.line = element_line(colour = 'black'),
            legend.position = "top",          # place legend on the right
            legend.title = element_blank(),     # optional: no title
            legend.text = element_text(size = 8,face = "bold"),
            axis.text = element_text(size = 12),
            axis.title = element_text(size = 14, face = "bold")) +
      # Apply custom color scheme
      scale_color_manual(values = custom_colors)+
      scale_fill_manual(values = custom_colors)
    
  } else {
    pcoa <- ggplot() +
      # Smaller scatter points
      geom_point(data = coordinates,
                 mapping = aes(x = `Axis.1`, y = `Axis.2`, colour = !!sym(grp)),
                 alpha = 0.6, size = 0.5) + 
      
      labs(title = paste0(dm_name, ': ', grp),
           subtitle = subtitle) +
      
      xlab(paste('PC1 (', label_percent()(pc_exp[1]), ')', sep = '')) +
      ylab(paste('PC2 (', label_percent()(pc_exp[2]), ')', sep = '')) +
      
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            axis.line = element_line(colour = 'black'),
            legend.position = "top",          # place legend on the right
            legend.title = element_blank(),     # optional: no title
            legend.text = element_text(size = 8,face = "bold"),
            axis.text = element_text(size = 12),
            axis.title = element_text(size = 14, face = "bold")) +
      # Apply custom color scheme
      scale_color_manual(values = custom_colors)+
      scale_fill_manual(values = custom_colors)
    
    
    
    #   pcoa <- ggplot() +
    #     geom_point(data = coordinates,
    #                mapping = aes(x = `Axis.1`, y = `Axis.2`, colour = !!sym(grp)),
    #                alpha = 0.8, size = 1) + 
    #     labs(title = paste0(dm_name, ': ', grp),
    #          subtitle = subtitle) +
    #     geom_segment(data = coordinates,
    #                  aes(x = CentroidX, xend = Axis.1, y = CentroidY, yend = Axis.2, color = !!sym(grp)),
    #                  alpha = 0.2, size = 0.5) +
    #     xlab(paste('PC1 (', label_percent()(pc_exp[1]), ')', sep = '')) +
    #     ylab(paste('PC2 (', label_percent()(pc_exp[2]), ')', sep = '')) +
    #     theme(panel.grid.major = element_blank(),
    #           panel.grid.minor = element_blank(),
    #           panel.background = element_blank(),
    #           axis.line = element_line(colour = 'black'),
    #           legend.title = element_blank(),
    #           legend.position = 'bottom',
    #           legend.text = element_text(size = 8)) +
    #     scale_color_gradient(low = "blue", high = "red")  # Adjust colors for your preference
    #   
  }
  
  print(pcoa)
  return(pcoa)
}

# function extracted from plot_unifrac to test adonis and permdisp
perm_test_customized_groups = function(dm,metatable,grp,dm_name){
    set.seed(666)
    # dummy input variable for testing function
    # dm = dm
    # metatable = metatable
    # grp = 'treatment'
    # 
    metatable = subset(metatable, get(grp) != 'N.A.', select = c(grp))
    dm = dist_subset(dm,rownames(metatable))
    stopifnot(all(rownames(metatable) == names(dm)))
    
    # Adonis test (only run if all groups have >= 5 samples)
    adonis.r = adonis2(dm~get(grp),metatable,permut = 1000,parallel = 1)
    #adonis.p = round(adonis.r$`Pr(>F)`[1],digit = 4)
    adonis.p = adonis.r$`Pr(>F)`[1]
    
    betadisp.r = permutest(betadisper(dm,metatable[[grp]]),permutations = 1000,parallel = 1)
    #betadisp.p = round(betadisp.r$tab[[6]][1],digit = 4)
    betadisp.p = betadisp.r$tab[[6]][1]
    
    cat('\n\n')
    print(capture.output(adonis.r))
    cat('\n\n')
    print(capture.output(betadisp.r))
    cat('\n\n')
    
    return(c(paste0(unique(metatable[[grp]]),collapse = ' VS '),adonis.p,betadisp.p))
}

# function for easy plotting pcoa plots (label each dot by patient ID)
plot_unifrac_labelbyID = function(dm,metatable,grp,dm_name,permtest = F,grp_continuous = F,custom_colors){
  # # dummy input variable for testing function
  # dm = dist.wUF
  # metatable = metatable
  # dm_name = 'x'
  # grp = 'disease_severity'
  
  metatable.full = metatable
  metatable = subset(metatable, get(grp) != 'N.A.', select = c(grp))
  # select distance matrix matching to metadata
  dm = dist_subset(dm,rownames(metatable))
  # some distance measures may result in negative eigenvalues. In that case, add a correction:
  pcoa.wunifrac = pcoa(dm,correction = 'cailliez')
  # transform the pcoa.wunifrac object into a df for making it plottable with ggplot2
  coordinates.pcoa.wunifrac = data.frame(pcoa.wunifrac$vectors,check.names = FALSE)
  
  # Adonis test (only run if all groups have >= 5 samples)
  ## Calculate the number of samples in each group
  group_counts <- table(metatable[[grp]])
  ## Check if any group has fewer than 5 samples
  any_small_groups <- any(group_counts < 5)
  if (any_small_groups == F & permtest ==T){
    stopifnot(all(rownames(metatable) == names(dm)))
    adonis.r = adonis2(dm~get(grp),metatable,permut = 1000,parallel = 1)
    adonis.p = round(adonis.r$`Pr(>F)`[1],digit = 4)
    
    betadisp.r = permutest(betadisper(dm,metatable[[grp]]),permutations = 1000,parallel = 1)
    betadisp.p = round(betadisp.r$tab[[6]][1],digit = 4)
    
    subtitle = paste0('permanova p val = ', adonis.p, ',betadisp p val = ', betadisp.p)
  } else {
    subtitle = 'Adonis/permdisp skipped'
  }
  
  
  # plot
  if (!all(rownames(metatable) == rownames(coordinates.pcoa.wunifrac))) {stop("mismatchs in pca coordinates and coloring column")} # sanity check
  
  # Calculate group centroids
  coordinates = coordinates.pcoa.wunifrac[,c(1,2)]
  coordinates$sample.id.new = rownames(coordinates)
  #coordinates[[grp]]= metatable[[grp]]      
  coordinates[[grp]] = metatable[[grp]][match(rownames(coordinates),rownames(metatable))] 
  #print(coordinates)
  centroids <- coordinates %>%
    group_by_at(grp) %>%
    summarize(CentroidX = mean(Axis.1), CentroidY = mean(Axis.2))
  coordinates = merge(coordinates,centroids,
                      by = grp)           
  coordinates$ID = metatable.full$ID[match(coordinates$sample.id.new,metatable.full$sample.id.new)]
  
  # automatically grab relative eigenvalue
  if ("Rel_corr_eig" %in% colnames(pcoa.wunifrac$values)) {
    pc_exp <- pcoa.wunifrac$values$Rel_corr_eig
  } else {
    pc_exp <- pcoa.wunifrac$values$Relative_eig
  }
  
  if (grp_continuous == F){
    
    pcoa <- ggplot() +
      # Smaller scatter points
      geom_point(data = coordinates,
                 mapping = aes(x = `Axis.1`, y = `Axis.2`, colour = get(grp)),
                 alpha = 0.6, size = 0.5) + 
      # Add text labels for each point
      geom_text(data = coordinates,
                aes(x = `Axis.1`, y = `Axis.2`, label = ID),
                size = 2, vjust = -0.5, check_overlap = TRUE) +
      # Connection lines to centroids
      geom_segment(data = coordinates,
                   aes(x = CentroidX, xend = Axis.1, y = CentroidY, yend = Axis.2, color = get(grp)),
                   alpha = 0.3, size = 0.3) +
      # Add ellipses
      #stat_ellipse(data = coordinates,
      #              aes(x = `Axis.1`, y = `Axis.2`, color = get(grp)),
      #               type = "norm", level = 0.68, size = 0.5, alpha = 0.4) +
      # Emphasize centroids
      geom_point(data = centroids,
                 aes(x = CentroidX, y = CentroidY, color = get(grp)),
                 size = 3, shape = 18, stroke = 1.2) +
      
      labs(title = paste0(dm_name, ': ', grp),
           subtitle = subtitle) +
      
      xlab(paste('PC1 (', label_percent()(pc_exp[1]), ')', sep = '')) +
      ylab(paste('PC2 (', label_percent()(pc_exp[2]), ')', sep = '')) +
      
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            axis.line = element_line(colour = 'black'),
            legend.position = "right",          # place legend on the right
            legend.direction = "vertical",      # make labels stack vertically
            legend.title = element_blank(),     # optional: no title
            legend.text = element_text(size = 8)) +
      # Apply custom color scheme
      scale_color_manual(values = custom_colors)
    
  } else {
    pcoa <- ggplot() +
      # Smaller scatter points
      geom_point(data = coordinates,
                 mapping = aes(x = `Axis.1`, y = `Axis.2`, colour = !!sym(grp)),
                 alpha = 0.6, size = 0.5) + 
      
      labs(title = paste0(dm_name, ': ', grp),
           subtitle = subtitle) +
      
      xlab(paste('PC1 (', label_percent()(pc_exp[1]), ')', sep = '')) +
      ylab(paste('PC2 (', label_percent()(pc_exp[2]), ')', sep = '')) +
      
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            axis.line = element_line(colour = 'black'),
            legend.position = "right",          # place legend on the right
            legend.direction = "vertical",      # make labels stack vertically
            legend.title = element_blank(),     # optional: no title
            legend.text = element_text(size = 8)) +
      # Apply custom color scheme
      scale_color_manual(values = custom_colors)  
    
    
    
    #   pcoa <- ggplot() +
    #     geom_point(data = coordinates,
    #                mapping = aes(x = `Axis.1`, y = `Axis.2`, colour = !!sym(grp)),
    #                alpha = 0.8, size = 1) + 
    #     labs(title = paste0(dm_name, ': ', grp),
    #          subtitle = subtitle) +
    #     geom_segment(data = coordinates,
    #                  aes(x = CentroidX, xend = Axis.1, y = CentroidY, yend = Axis.2, color = !!sym(grp)),
    #                  alpha = 0.2, size = 0.5) +
    #     xlab(paste('PC1 (', label_percent()(pc_exp[1]), ')', sep = '')) +
    #     ylab(paste('PC2 (', label_percent()(pc_exp[2]), ')', sep = '')) +
    #     theme(panel.grid.major = element_blank(),
    #           panel.grid.minor = element_blank(),
    #           panel.background = element_blank(),
    #           axis.line = element_line(colour = 'black'),
    #           legend.title = element_blank(),
    #           legend.position = 'bottom',
    #           legend.text = element_text(size = 8)) +
    #     scale_color_gradient(low = "blue", high = "red")  # Adjust colors for your preference
    #   
  }
  
  print(pcoa)
}


# Plot ANCOMBC results
plot_ancombc = function(pu.raw,grp) {
  # remove N.A. samples from the desired grouping column
  pu = subset_samples(pu.raw,get(grp) != 'N.A.')
  # iterate the same process to test on genus and ASV level
  for (type in c('Genus','ASV')){
    # run ancombc
    if (type == 'Genus'){
      # the only difference between the two runs is that if running on ASV level, argument tax_level = NULL
      pu.ac= ancombc(phyloseq = pu, #the input data (phyloseq object)
                     formula = grp, #the character string expresses how microbial absolute abundances for each taxon depend on the variables in metadata.
                     tax_level = "Genus", #character. The taxonomic level of interest. The input data can be agglomerated at different taxonomic levels
                     p_adj_method = "holm", #character. method to adjust p-values. Default is "holm".
                     prv_cut = 0.1, #a numerical fraction between 0 and 1. Taxa with prevalences less than prv_cut will be excluded in the analysis. Default is 0.10.
                     lib_cut = 0, #a numerical threshold for filtering samples based on library sizes. Samples with library sizes less than lib_cut will be excluded in the analysis. Default is 0, i.e. do not discard any sample.
                     group = NULL, #character. the name of the group variable in metadata. group should be discrete. Specifying group is required for detecting structural zeros and performing global test. Default is NULL. If the group of interest contains only two categories, leave it as NULL.
                     struc_zero = F, #logical. whether to detect structural zeros based on group. Default is FALSE.
                     neg_lb = F, #logical. whether to classify a taxon as a structural zero using its asymptotic lower bound. Default is FALSE.
                     tol = 1e-5, #numeric. the iteration convergence tolerance for the E-M algorithm. Default is 1e-05.
                     max_iter = 100, #numeric. the maximum number of iterations for the E-M algorithm. Default is 100.
                     conserve = FALSE, #logical. whether to use a conservative variance estimator for the test statistic. It is recommended if the sample size is small and/or the number of differentially abundant taxa is believed to be large. Default is FALSE.
                     alpha = 0.05, #numeric. level of significance. Default is 0.05.
                     global = F, #logical. whether to perform the global test. Default is FALSE.
                     n_cl = 1, #numeric. The number of nodes to be forked. For details, see ?parallel::makeCluster. Default is 1 (no parallel computing).
                     verbose = F) #logical. Whether to generate verbose output during the ANCOM-BC fitting process. Default is FALSE.
    } else {
      pu.ac= ancombc(phyloseq = pu, #the input data (phyloseq object)
                     formula = grp, #the character string expresses how microbial absolute abundances for each taxon depend on the variables in metadata.
                     tax_level = NULL, #character. The taxonomic level of interest. The input data can be agglomerated at different taxonomic levels
                     p_adj_method = "holm", #character. method to adjust p-values. Default is "holm".
                     prv_cut = 0.1, #a numerical fraction between 0 and 1. Taxa with prevalences less than prv_cut will be excluded in the analysis. Default is 0.10.
                     lib_cut = 0, #a numerical threshold for filtering samples based on library sizes. Samples with library sizes less than lib_cut will be excluded in the analysis. Default is 0, i.e. do not discard any sample.
                     group = NULL, #character. the name of the group variable in metadata. group should be discrete. Specifying group is required for detecting structural zeros and performing global test. Default is NULL. If the group of interest contains only two categories, leave it as NULL.
                     struc_zero = F, #logical. whether to detect structural zeros based on group. Default is FALSE.
                     neg_lb = F, #logical. whether to classify a taxon as a structural zero using its asymptotic lower bound. Default is FALSE.
                     tol = 1e-5, #numeric. the iteration convergence tolerance for the E-M algorithm. Default is 1e-05.
                     max_iter = 100, #numeric. the maximum number of iterations for the E-M algorithm. Default is 100.
                     conserve = FALSE, #logical. whether to use a conservative variance estimator for the test statistic. It is recommended if the sample size is small and/or the number of differentially abundant taxa is believed to be large. Default is FALSE.
                     alpha = 0.05, #numeric. level of significance. Default is 0.05.
                     global = F, #logical. whether to perform the global test. Default is FALSE.
                     n_cl = 1, #numeric. The number of nodes to be forked. For details, see ?parallel::makeCluster. Default is 1 (no parallel computing).
                     verbose = F) #logical. Whether to generate verbose output during the ANCOM-BC fitting process. Default is FALSE.
    }
    
    # process result table
    res = as.data.frame(pu.ac$res)
    # remove unwanted columns  
    res = res %>%
      select(-c(se.taxon,W.taxon,p_val.taxon,q_val.taxon,diff_abn.taxon, 
                names(res)[grepl('Intercept',names(res))]))
    names(res)[names(res) == 'lfc.taxon'] = 'taxon'
    
    # stack it (even if there are only 2 groups compared, which results in no columns to stack, this line still works and functions as the way to separate column names based on ".grp", and take whatever after the separator to the new column "Test_level")
    res <- pivot_longer(res, cols = -taxon, names_to = c(".value", "Test_level"), names_sep = paste0(".", grp))
    
    # if it is genus level comparison
    if (type == 'Genus'){
      res = subset(res, grepl('Genus',taxon))
      res$taxon = gsub('Genus:','',res$taxon)
      res = merge(res, distinct(select(as.data.frame(tax_table(pu)),c(Genus,Family,Phylum))),
                  by.x = 'taxon', by.y = 'Genus',
                  all.x = T)
      
    } else {
      res = merge(res, select(as.data.frame(tax_table(pu)),c(Family,Phylum)),
                  by.x = 'taxon', by.y = 0,
                  all.x = T)
      # Some ASVs don't have taxonomy classified, let's remove them at this step
      res = res[complete.cases(res),]
    }
    
    
    res.p.05 = subset(res,p_val < 0.05)
    
    # get the reference level in ANCOMBC testing
    ref_level = setdiff(unique(sample_data(pu)[[grp]]), unique(res$Test_level))
    p = ggplot(res.p.05,aes(x=Phylum,y=lfc,color=Family,shape = factor(diff_abn)))+
      geom_jitter(size=3,alpha=0.7,width=0.1) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
      labs(title = paste0('ANCOMBC at ', type, ' level : testing ', grp, " with ", ref_level, ' as reference'),
           subtitle = paste0("Only showing ", type, " with p < 0.05"),
           color = "Family",
           shape = 'Test significance (q < 0.05)') +
      labs(x = "Phylum", y = "Log Fold Change")
    # add facet if there are more than 1 comparison (more than 2 groups compared)
    if (length(unique(res.p.05$Test_level))>1) {
      p = p + facet_wrap(~Test_level)
    }
    print(p)
    assign(paste0('plt.',type),p)
    assign(paste0('res.',type),res)
  }
  return(list(res.ASV,res.Genus,plt.ASV,plt.Genus))
}

# count total number of reads for each sample and add it as a column to phyloseq metadata
count_total_reads = function(psq,suffix){
  total_reads = as.data.frame(sample_sums(psq))
  names(total_reads) = 'total_reads'
  sample_data(psq)[[paste0('read_count_',suffix)]] = total_reads$total_reads[match(rownames(sample_data(psq)),rownames(total_reads))]
  return(psq)
}

create_dt <- function(df, caption) {
  datatable(df,
            caption = caption,
            extensions = 'Buttons',
            options = list(dom = 'Blfrtip',
                           buttons = c('copy', 'csv', 'excel', 'pdf', 'print'),
                           lengthMenu = list(c(10, 25, 50, -1),
                                             c(10, 25, 50, "All")),
                           pageLength = 20, # Set the default rows per page here
                           scrollY = 400,    # Set vertical scroll height to 400 pixels (adjust as needed)))  
                           scrollCollapse = TRUE))
}

# function to count patients with samples taking at different pairs of stages; taking a metadata dataframe  as the input
# depth argument is used to customize the column showing the sample count for easy merging afterwards
calculate_stage_combinations <- function(df,depth) {
  # df <- psq %>%
  #   sample_data %>%
  #   as.matrix %>%
  #   as.data.frame
  
  # Check if columns ID and Stage exist in the data frame
  if(!all(c("ID", "Stage") %in% colnames(df))) {
    stop("Data frame must contain columns named 'ID' and 'Stage'")
  }
  
  # Aggregate stages for each ID, creating a unique combination of stages per ID
  id_combinations <- aggregate(Stage ~ ID, data = df, FUN = function(x) paste(sort(unique(x)), collapse = ", "))
  
  # Count the number of IDs with each unique stage combination
  combination_counts <- table(id_combinations$Stage)
  
  # Convert to a data frame for easier readability
  combination_counts_df <- as.data.frame(combination_counts)
  colnames(combination_counts_df) <- c("Stage_Combination", "Count")
  
  # Count the possible paired comparisons
  combination_counts_df$Stage_Combination = as.character(combination_counts_df$Stage_Combination)
  # Define paired comparisons
  paired_comparisons <- list(
    '1_paired_comparison_A_R' = c("A", "R"),
    '2_paired_comparison_R_C' = c("R", "C"),
    '3_paired_comparison_C_C24' = c("C", "C24"),
    '4_paired_comparison_A_C' = c("A", "C"),
    '5_paired_comparison_A_C24' = c("A", "C24"),
    '6_paired_comparison_C24_R' = c("C24", "R")
  )
  
  # Function to calculate the sum of counts for each paired comparison
  get_sum <- function(keywords, df) {
    sum(df$Count[sapply(df$Stage_Combination, function(x) all(keywords %in% unlist(strsplit(x, ",\\s*"))))])
  }
  
  # Calculate sums for each paired comparison
  results <- sapply(paired_comparisons, get_sum, df = combination_counts_df)
  
  # Create the new DataFrame
  new_df <- data.frame(Stage_Combination = names(results), Count = results)
  # merge it with the previous one
  combination_counts_df = bind_rows(new_df,combination_counts_df)
  # customize the name of the count column for easy merging
  names(combination_counts_df)[names(combination_counts_df) == 'Count'] <- paste0('RRF_at_',depth)
  return(combination_counts_df)
}

# Agglomerate a phyloseq object to a specific taxrank.
# input is a phyloseq obj, out put a dataframe with rows as samples and columns as taxrank selected
tax_agg = function(p,sampleindex = 'Sample', taxrank = 'Genus',remove_unclassified = F){
  tg <- p %>%
    psmelt %>%
    group_by_at(c(sampleindex,taxrank)) %>%
    summarise(x = sum(Abundance)) %>%
    pivot_wider(names_from = taxrank, values_from = x) %>%
    as.data.frame %>%
    setrowname
  if (remove_unclassified) {
    tg = select(tg,-Unclassified)
  }
  return(tg)
}

# function for easy plotting PCoA + envfit
# dm metrics and envfitDf do not need to be subset. envfitDf needs to have rows as samples with rownames exactly match to dm and metatable rowname

plot_envifit = function(dm,metatable,grp,dm_name,envfitDf,topN_envfit_variables_to_plot,
                        permtest = F,arrowshrinkingby=0.5,custom_colors,plot_significant_arrow_only = T){
  # # dummy input variable for testing function
  # grp = 'Stage'
  # dm = dist.wUF
  # metatable <- p.all.cpm %>%
  #   sample_data %>%
  #   as.matrix %>%
  #   as.data.frame
  # dm_name = 'Weighted UniFrac'
  # envfitDf2 = tax_agg(p.all.cpm,remove_unclassified = F)
  # topN_envfit_variables_to_plot = 5
  # permtest = T
  # xnudge=0.2
  # ynudge=0.2
  # 
  set.seed(666)
  metatable = subset(metatable, get(grp) != 'N.A.', select = c(grp))
  # select distance matrix matching to metadata
  dm = dist_subset(dm,rownames(metatable))
  # some distance measures may result in negative eigenvalues. In that case, add a correction:
  pcoa.wunifrac = pcoa(dm,correction = 'cailliez')
  # transform the pcoa.wunifrac object into a df for making it plottable with ggplot2
  coordinates.pcoa.wunifrac = data.frame(pcoa.wunifrac$vectors,check.names = FALSE)
  
  # Adonis test (only run if all groups have >= 5 samples)
  ## Calculate the number of samples in each group
  group_counts <- table(metatable[[grp]])
  ## Check if any group has fewer than 5 samples
  any_small_groups <- any(group_counts < 5)
  if (any_small_groups == F & permtest == T){
    stopifnot(all(rownames(metatable) == names(dm)))
    adonis.r = adonis2(dm~get(grp),metatable,permut = 1000,parallel = 1)
    adonis.p = round(adonis.r$`Pr(>F)`[1],digit = 4)
    
    betadisp.r = permutest(betadisper(dm,metatable[[grp]]),permutations = 1000,parallel = 1)
    betadisp.p = round(betadisp.r$tab[[6]][1],digit = 4)
    
    subtitle = paste0('p(adonis) = ', adonis.p, ', p(betadisp) = ', betadisp.p)
  } else {
    subtitle = 'Adonis/permdisp skipped'
  }
  
  
  # plot
  if (!all(rownames(metatable) == rownames(coordinates.pcoa.wunifrac))) {stop("mismatchs in pca coordinates and coloring column")} # sanity check
  
  # Calculate group centroids
  coordinates = coordinates.pcoa.wunifrac[,c(1,2)]
  #coordinates[[grp]]= metatable[[grp]]
  coordinates[[grp]] = metatable[[grp]][match(rownames(coordinates),rownames(metatable))] 
  coordinates$sampleid = rownames(coordinates)
  centroids <- coordinates %>%
    group_by_at(grp) %>%
    summarize(CentroidX = mean(Axis.1), CentroidY = mean(Axis.2))
  coordinates = merge(coordinates,centroids,
                      by = grp)           
  rownames(coordinates) = coordinates$sampleid
  coordinates = coordinates[order(coordinates$sampleid),]
  
  # automatically grab relative eigenvalue
  if ("Rel_corr_eig" %in% colnames(pcoa.wunifrac$values)) {
    pc_exp <- pcoa.wunifrac$values$Rel_corr_eig
  } else {
    pc_exp <- pcoa.wunifrac$values$Relative_eig
  }
  
  # envfit calculation
  # match samples between pcoa coordinates and envfit table
  envfitDf = envfitDf[rownames(coordinates),]
  # remove the unclassified ones
  envfitDf <- envfitDf %>%
    dplyr::select(-matches("unclassified"))
  stopifnot(rownames(coordinates) == rownames(envfitDf))
  # perform standardization column-wisely (each column is a variable/feature/taxon)
  envfitDf = decostand(envfitDf, method = "standardize", MARGIN = 2)
  # sometimes scaling generates NA, especially when data is very sparse; so lets remove it
  envfitDf <- envfitDf[, colSums(is.na(envfitDf)) == 0]
  
  envfit_result <- envfit(dplyr::select(coordinates,c(Axis.1,Axis.2)), envfitDf, perm = 999)
  # Fortify the envfit results 
  vectors <- as.data.frame(scores(envfit_result, "vectors"))  # Extract vector coordinate
  vectors$r2 <- envfit_result$vectors$r  # Add R² values
  vectors$pval <- envfit_result$vectors$pvals  # Add p-values
  
  # FDR correction
  ## Calculate FDR using Benjamini-Hochberg method and add it as a new column
  vectors$FDR <- p.adjust(vectors$pval, method = "BH")
  
  # Add 'sig' column based on FDR value ranges
  vectors <- vectors %>%
    mutate(sig = case_when(
      FDR < 0.001 ~ "FDR < 0.001",
      FDR >= 0.001 & FDR < 0.01 ~ "FDR < 0.01",
      FDR >= 0.01 & FDR < 0.05 ~ "FDR < 0.05",
      FDR >= 0.05 ~ "ns"
    ))
  
  # Select top x vectors with the highest R2 
  vectors_s <- vectors %>%
    arrange(desc(r2)) %>%
    slice_head(n = topN_envfit_variables_to_plot)
  
  if (plot_significant_arrow_only) {
    vectors_s <- vectors_s %>% filter(FDR < 0.05)
  }
  # Shrink arrow size to fit the PCoA plot
  vectors_s$Axis.1 = vectors_s$Axis.1 * arrowshrinkingby
  vectors_s$Axis.2 = vectors_s$Axis.2 * arrowshrinkingby
  # Calculate x and y axis limits based on the range of data points and expand by 10%
  x_range <- range(c(coordinates$Axis.1, vectors_s$Axis.1), na.rm = TRUE)
  y_range <- range(c(coordinates$Axis.2, vectors_s$Axis.2), na.rm = TRUE)
  
  # Expand limits by 10%
  x_limits <- c(x_range[1] , x_range[2])
  y_limits <- c(y_range[1] , y_range[2] )
  
  pcoa <- ggplot() +
    # Smaller scatter points
    geom_point(data = coordinates,
               mapping = aes(x = `Axis.1`, y = `Axis.2`, color = get(grp)),
               alpha = 0.6, size = 0.5) + 
    
    # Connection lines to centroids
    geom_segment(data = coordinates,
                 aes(x = CentroidX, xend = Axis.1, y = CentroidY, yend = Axis.2, color = get(grp)),
                 alpha = 0.3, size = 0.3) +
    # Add ellipses
    #stat_ellipse(data = coordinates,
    #              aes(x = `Axis.1`, y = `Axis.2`, color = get(grp)),
    #               type = "norm", level = 0.68, size = 0.5, alpha = 0.4) +
    # Emphasize centroids
    geom_point(data = centroids,
               aes(x = CentroidX, y = CentroidY, fill = get(grp)),
               size = 4, shape = 21, stroke = 0.6, color = 'black') +
    # envifit
    geom_segment(data = vectors_s,
                 aes(x = 0, y = 0, xend = Axis.1, yend = Axis.2, linetype = sig),
                 color = 'blue',alpha = 1.2,
                 arrow = arrow(length = unit(0.3, "cm")), size = 0.5, show.legend = FALSE) +
    
    # Axis limits
    xlim(x_limits) +
    ylim(y_limits) +
    
    xlab(paste('PC1 (',label_percent()(pc_exp[1]),')',sep = '')) +
    ylab(paste('PC2 (',label_percent()(pc_exp[2]),')',sep = '')) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = 'black'),
          legend.position = "top",          
          legend.title = element_blank(),     # optional: no title
          legend.text = element_text(size = 8,face = "bold"),
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 14, face = "bold")) +
    scale_fill_manual(values = custom_colors) +
    scale_color_manual(values = custom_colors) +
    guides(fill= guide_legend(nrow=2))  
    # scale_linetype_manual(name = "Significance (FDR)",
    #                       values = c("FDR < 0.001" = "dotdash", 
    #                                  "FDR < 0.01" = "dashed", 
    #                                  "FDR < 0.05" = "solid", 
    #                                  "ns" = "dotted"))
  pcoa_ggtextrepel= pcoa + geom_text_repel(data = vectors_s, aes(x = Axis.1, y = Axis.2, label = rownames(vectors_s),
                                                                 nudge_x = ifelse(Axis.1 > 0, 0.3, -0.3),  # Nudge right if x > 0, left if x < 0
                                                                 nudge_y = ifelse(Axis.2 > 0, 0.3, -0.3)), # Nudge up if y > 0, down if y < 0
                                           box.padding = 0.6,
                                           point.padding = 0.4,
                                           segment.curvature = -0.1,
                                           segment.ncp = 3,
                                           force = 6,            # <- push harder
                                           force_pull = 0.05,    # <- pull back less
                                           max.time = 2,         # <- allow more optimization time (if available)
                                           # segment.angle = 10,
                                           direction = "both",
                                           max.overlaps = Inf,
                                           colour = 'blue', alpha = 0.9
  ) 
  pcoa_ggtext = pcoa + geom_text(data = vectors_s, aes(x = Axis.1, y = Axis.2, label = rownames(vectors_s),colour = 'blue'),
                                 nudge_x = 0,         # Slightly shift text horizontally
                                 nudge_y = 0,         # Slightly shift text vertically
                                 check_overlap = TRUE    # Prevent overlapping labels
  )
  return(list(pcoa_ggtextrepel,pcoa_ggtext,coordinates,vectors))
}

# execute permanova and permdisp specifically designed for testing Stage pairs
perm_test = function(metatable,dm,groups,sectionLevel = '#### '){
  set.seed(666)
  # select samples
  metatable.s = subset(metatable,Stage %in% groups)
  dm = dist_subset(dm,rownames(metatable.s))
  stopifnot(rownames(as.matrix(dm)) == rownames(metatable.s))

  # adonis
  adonis.r = adonis2(dm~Stage,metatable.s,permut = 1000, parallel = 1)
  adonis.p = adonis.r$`Pr(>F)`[1]
  # betadisp
  betadisp.r = permutest(betadisper(dm,metatable.s[['Stage']]),permutations = 1000,parallel = 1)
  betadisp.p = betadisp.r$tab[[6]][1]
  
  cat('\n\n')
  cat(paste0(sectionLevel,paste0(groups,collapse = ' VS ')))
  cat('\n\n')
  print(capture.output(adonis.r))
  cat('\n\n')
  print(capture.output(betadisp.r))
  cat('\n\n')
  
  return(c(paste0(groups,collapse = ' VS '),adonis.p,betadisp.p))
}

# execute permanova and permdisp specifically designed for testing stage pairs but block by patient
perm_test.blocked = function(metatable,dm,groups,sectionLevel = '#### '){
  set.seed(666)
  # dummy variables
  # dm = dist.wUF
  # metatable <- p.1r %>%
  #   sample_data %>%
  #   as.matrix %>%
  #   as.data.frame
  # groups =c('Convalescent','Recovery')
  
  # select samples (select the paired samples)
  metatable$rowname = rownames(metatable)
  # print(metatable)
  metatable.s <- metatable %>%
    subset(Stage %in% groups) %>%
    group_by(ID) %>%
    filter(n_distinct(Stage) == 2) %>%
    ungroup() %>%
    # Arrange by Patient.ID and Stage to ensure correct pairing order
    arrange(ID, Stage)
  rownames(metatable.s) = metatable.s$rowname
  
  # count the number of patients in the permutation
  patient.count = length(unique(metatable.s$ID))
  
  # subset dm
  dm = dist_subset(dm,rownames(metatable.s))
  stopifnot(rownames(as.matrix(dm)) == rownames(metatable.s))
  
  # block on patient ID
  h1 <- with(metatable.s, how(nperm = 999, blocks = ID))
  
  # adonis
  adonis.r = adonis2(dm~Stage, data = metatable.s, permutations = h1, parallel = 1)
  adonis.p = adonis.r$`Pr(>F)`[1]
  # betadisp
  betadisp.r = permutest(betadisper(dm,metatable.s[['Stage']]),permutations = h1, parallel = 1)
  betadisp.p = betadisp.r$tab[[6]][1]
  
  cat('\n\n')
  cat(paste0(sectionLevel,paste0(groups,collapse = ' VS ')))
  cat('\n\n')
  print(capture.output(adonis.r))
  cat('\n\n')
  print(capture.output(betadisp.r))
  cat('\n\n')
  
  return(c(paste0(paste0(groups,collapse = ' VS '), ' blocked by PatientID'),adonis.p,betadisp.p,patient.count))
}

# use phyloseq function merge_samples to merge replicated samples
## The function first extract samples based on sample index given as a argument from a phylosesq obj. The extracted phyloseq object will be merged
## (please note that it assumes all samples within the extracted one will be merged, so when customizing metadata values, it follows this rational)
## It then merge samples by 1) OTU table will be a sum of the given sample (which will be checked in the end, assuming only two samples are merged).
## 2) sample_data will be 1) if character column, return NA (in this function, designated column will be assigned value as shown below) 2) if numeric column, it will calculate average by default
## The function returns a merged phyloseq object

merge_duplicates = function(sampleIndexVector,obj = p.all){
  
  psq = prune_samples(sampleIndexVector,obj)
  y = as.data.frame(otu_table(psq))
  
  psq2 = merge_samples(psq,'ID_Stage')
  sample_data(psq2)[,c('uniqID','sample.id')] = c(paste0(sample_data(psq)[['uniqID']], collapse = '_and_'),paste0(sample_data(psq)[['sample.id']],collapse = '_and_'))
  sample_data(psq2)[,c('specimen_control','specimen_type','ID','Date','Stage','ID_Stage')] = 
    sample_data(psq)[1,c('specimen_control','specimen_type','ID','Date','Stage','ID_Stage')]
  sample_data(psq2)[,c('read_count_raw','read_count_before_rarefaction')] = 
    c(sum(sample_data(psq)[['read_count_raw']]),sum(sample_data(psq)[['read_count_before_rarefaction']]))
  y2 = as.data.frame(t(as.matrix(otu_table(psq2))))
  y3 = merge(y,y2,by=0)
  stopifnot(y3[,2] + y3[,3] == y3[,4])
  
  x3 = as.data.frame(t(as.matrix(bind_rows(as.data.frame(as.matrix(sample_data(psq))),as.data.frame(as.matrix(sample_data(psq2)))))))
  # print(x3)
  # print(y3)
  return(psq2)
}

# get top 10 taxa at the given taxonomic level from a phyloseq object
taxReadCount_collapse_top10_from_phyloseq = function(psq,tax,sample_list){
  # 
  # psq <- psq.tmp
  # tax <- 'Species'
  # sample_list <- sample_names(psq.tmp)
  intable = merge(as.data.frame(as.matrix(tax_table(psq))),
                  as.data.frame(as.matrix(otu_table(psq))),
                  by = 0)
  outtable = NULL
  outtable.10 = NULL
  outtable.10s = NULL
  # output is outtable
  outtable =select(intable, c(tax,sample_list))
  outtable = outtable %>%
    group_by_at(vars(tax)) %>%
    summarise_all(sum)
  # order the rows from high to low
  outtable = outtable[order(rowSums(outtable[,-1]), decreasing = T),]
  
  if (all(round(colSums(outtable[,-1])) == round(colSums(select(intable, c(sample_list)))))){
    # collapse the table to 9 most abundant entries and "the others"
    # output table is outtable.10
    if (nrow(outtable) >10) {
      # bind 1-9 to the sum of all others
      outtable.10 = bind_rows(outtable[1:9,],colSums(outtable[10:nrow(outtable),-1]))
      outtable.10[10,1] = 'Others'
    } else (outtable.10 = outtable)
    
    # stack the outtable.10
    # output is outtable.10s
    if (all(round(colSums(outtable.10[,-1]))== round(colSums(outtable[,-1])))){
      # stack the table
      outtable.10s = gather(outtable.10,"sample",'frequency',-tax)
      # # create a grouping variable
      # outtable.10s$group = ifelse(startsWith(outtable.10s$sample,'Control'),'Control',study)
      # outtable.10s$group = as.factor(outtable.10s$group)
    } else {outtable.10 <- outtable.10s <- NULL}
    
    # sanity check for the stacked table
    if (round(sum(outtable.10s$frequency)) != round(sum(select(intable, c(sample_list))))){
      outtable.10s = NULL
    }
  } else {outtable <- outtable.10 <- outtable.10s <-  NULL}
  outtable_all = list(outtable, outtable.10, outtable.10s)
  return(outtable_all)
}

# get top 20 taxa at the given taxonomic level from a phyloseq object
taxReadCount_collapse_top20_from_phyloseq = function(psq,tax,sample_list){
  # 
  # psq <- psq.tmp
  # tax <- 'Genus'
  # sample_list <- sample_names(psq.tmp)
  intable = merge(as.data.frame(as.matrix(tax_table(psq))),
                  as.data.frame(as.matrix(otu_table(psq))),
                  by = 0)
  outtable = NULL
  outtable.10 = NULL
  outtable.10s = NULL
  # output is outtable
  outtable =select(intable, c(tax,sample_list))
  outtable = outtable %>%
    group_by_at(vars(tax)) %>%
    summarise_all(sum)
  # order the rows from high to low
  outtable = outtable[order(rowSums(outtable[,-1]), decreasing = T),]
  
  if (all(round(colSums(outtable[,-1])) == round(colSums(select(intable, c(sample_list)))))){
    # collapse the table to 19 most abundant entries and "the others"
    # output table is outtable.10
    if (nrow(outtable) >20) {
      # bind 1-9 to the sum of all others
      outtable.10 = bind_rows(outtable[1:19,],colSums(outtable[20:nrow(outtable),-1]))
      outtable.10[20,1] = 'Others'
    } else (outtable.10 = outtable)
    
    # stack the outtable.10
    # output is outtable.10s
    if (all(round(colSums(outtable.10[,-1]))== round(colSums(outtable[,-1])))){
      # stack the table
      outtable.10s = gather(outtable.10,"sample",'frequency',-tax)
      # # create a grouping variable
      # outtable.10s$group = ifelse(startsWith(outtable.10s$sample,'Control'),'Control',study)
      # outtable.10s$group = as.factor(outtable.10s$group)
    } else {outtable.10 <- outtable.10s <- NULL}
    
    # sanity check for the stacked table
    if (round(sum(outtable.10s$frequency)) != round(sum(select(intable, c(sample_list))))){
      outtable.10s = NULL
    }
  } else {outtable <- outtable.10 <- outtable.10s <-  NULL}
  outtable_all = list(outtable, outtable.10, outtable.10s)
  return(outtable_all)
}

# get top x taxa at the given taxonomic level from a phyloseq object
# if merge_unclassified = T, then all unclassified ones (like xxx_unclassified and yyy_unclassified) will be grouped into unclassified
# and when calculate top X taxa, the unclassified will be included into Others
taxReadCount_collapse_topn_from_phyloseq = function(psq,tax,sample_list,topn = 10, merge_unclassified = F){
  # 
  # psq <- psq.tmp
  # xtax <- 'Genus'
  # sample_list <- sample_names(psq.tmp)
  # topn = 10
  # merge_unclassified = T
  
  intable = merge(as.data.frame(as.matrix(tax_table(psq))),
                  as.data.frame(as.matrix(otu_table(psq))),
                  by = 0)
  outtable = NULL
  outtable.10 = NULL
  outtable.10s = NULL
  
  # output is outtable
  outtable =select(intable, c(tax,sample_list))
  
  if (merge_unclassified) {
    outtable[[tax]] = ifelse(grepl('unclassified',outtable[[tax]],ignore.case = T),'unclassified',outtable[[tax]])
  }
  
  outtable = outtable %>%
    group_by_at(vars(tax)) %>%
    summarise_all(sum)
  
  # order the rows from high to low
  # if ignore unclassified
  if (merge_unclassified){
    outtable <- outtable[
      order(outtable[[1]] == "unclassified",   # 1st sorting key
            -rowSums(outtable[,-1])            # 2nd sorting key
      ),
    ]
  } else {
    outtable = outtable[order(rowSums(outtable[,-1]), decreasing = T),]
  }
  
 
  
  if (all(round(colSums(outtable[,-1])) == round(colSums(select(intable, c(sample_list)))))){
    # collapse the table to x most abundant entries and "the others"
    # output table is outtable.10
    if (nrow(outtable) > topn + 1) {
      # bind 1-topn to the sum of all others
      outtable.10 = bind_rows(outtable[1:topn,],colSums(outtable[topn + 1:nrow(outtable),-1],na.rm = T))
      outtable.10[topn+1,1] = 'Others'
    } else (outtable.10 = outtable)
    
    # stack the outtable.10
    # output is outtable.10s
    if (all(round(colSums(outtable.10[,-1]))== round(colSums(outtable[,-1])))){
      # stack the table
      outtable.10s = gather(outtable.10,"sample",'frequency',-tax)
      # # create a grouping variable
      # outtable.10s$group = ifelse(startsWith(outtable.10s$sample,'Control'),'Control',study)
      # outtable.10s$group = as.factor(outtable.10s$group)
    } else {outtable.10 <- outtable.10s <- NULL}
    
    # sanity check for the stacked table
    if (round(sum(outtable.10s$frequency)) != round(sum(select(intable, c(sample_list))))){
      outtable.10s = NULL
    }
  } else {outtable <- outtable.10 <- outtable.10s <-  NULL}
  outtable_all = list(outtable, outtable.10, outtable.10s)
  return(outtable_all)
}


# take a preranked taxa vector, place "unclassified" and "others" at the end
taxa_RA_rerank = function(x){
  if ("Others" %in% x){
    x = c(setdiff(x,c('Others')),'Others')
  }
  if ('Unclassified' %in% x){
    x = c(setdiff(x,c('Unclassified')),'Unclassified')
  }
  return(x)
}

# wrapper to run maaslin2 with a phyloseq object. If feature is any taxonomy level (e.g. Genus), it automatically agglomerate to that level
# If it is OTU, it will use the OTU table as input.
# the other part of the function is hard coded for this project -- it will automatically extract ID and daySinceSymptomOnset from metadata,
# and run maaslin2 using ID as random effect and daySinceSymptomOnset as fixed effect
run_maaslin = function(psq.tmp,res_dir,feature){
  # psq.tmp = psq.tmp
  # #res_dir = paste0('../data/Maaslin2/daySinceSymptomOnset')
  # feature = 'Genus'
  
  # extract relative abundance matrix
  if (feature != 'OTU'){
    exps <- taxReadCount_collapse_top10_from_phyloseq(psq.tmp,feature,sample_list = sample_names(psq.tmp))[[1]] %>%
      as.matrix %>%
      as.data.frame %>%
      setrowname
  } else {
    exps <- otu_table(psq.tmp) %>%
      as.matrix %>%
      as.data.frame
  }
  exps[] <- lapply(exps, function(x) as.numeric(as.character(x)))
  if (!all(round(colSums(exps),digits = 1)==100)) {
    stop("Error: sample-wise sum of relative abundance != 100 for some samples")
  }
  # removed the unclassified features
  exps = subset(exps, !grepl('unclassified',rownames(exps)))
  
  # extract phenotype table
  pheno <- psq.tmp %>%
    sample_data %>%
    as.matrix %>%
    as.data.frame %>%
    select(c(ID,daysSinceSymptomOnset)) %>%
    mutate(daysSinceSymptomOnset = as.numeric(daysSinceSymptomOnset),
           ID = paste0('patient',ID))

  # sanity check
  if (!all(rownames(pheno) == colnames(exps))) {
    stop("Error: rownames of 'pheno' do not match colnames of 'exps'. Please check your data.")
  }

  # make a new directory for placing the result folder
  res_dir = res_dir
  if (file.exists(res_dir)){
    unlink(res_dir,recursive = TRUE)
    dir.create(res_dir,recursive = T)
  } else {
    dir.create(res_dir,recursive = T)
  }
  
  fit1 = capture.output(# capture.out is used to mutate verbose output by maaslin
    Maaslin2(input_data = exps, 
             input_metadata = pheno,
             min_prevalence = 0.5, # the minimum percent of samples for which a feature is detected at minimum abundance (default 0.1)
             min_abundance = 0.05, # The minimum abundance for each feature (default 0)
             normalization = 'None', # default TSS
             transform = 'LOG', # default log 
             fixed_effects = c('daysSinceSymptomOnset'),
             random_effects = "ID",
             cores = 8,
             max_significance = 0.3, # this is set to high value so that maaslin automatically plot every figure that might need in future, #FDR q value threshold
             output = res_dir))
  res = read.csv(paste(res_dir, 'all_results.tsv', sep = '/'), sep = '\t')
  return(list(exps,pheno,res))
}

run_maaslin_unrarefied = function(psq.tmp = psq.tmp,res_dir = res_dir,feature = feature){
  # psq.tmp = psq.tmp
  # res_dir = paste0('../data/Maaslin2/daySinceSymptomOnset')
  # feature = 'Genus'
  
  # extract relative abundance matrix
  if (feature != 'OTU'){
    exps <- taxReadCount_collapse_top10_from_phyloseq(psq.tmp,feature,sample_list = sample_names(psq.tmp))[[1]] %>%
      as.matrix %>%
      as.data.frame %>%
      setrowname
  } else {
    exps <- otu_table(psq.tmp) %>%
      as.matrix %>%
      as.data.frame
  }
  exps[] <- lapply(exps, function(x) as.numeric(as.character(x)))
  # if (!all(round(colSums(exps),digits = 1)==100)) {
  #   stop("Error: sample-wise sum of relative abundance != 100 for some samples")
  # }
  # removed the unclassified features
  exps = subset(exps, !grepl('unclassified',rownames(exps)))
  
  # extract phenotype table
  pheno <- psq.tmp %>%
    sample_data %>%
    as.matrix %>%
    as.data.frame %>%
    select(c(ID,daysSinceSymptomOnset))
  
  # pheno$gender = as.factor(pheno$gender)
  # pheno$race = as.factor(pheno$race)
  # pheno$age = as.numeric(pheno$age)
  pheno$daysSinceSymptomOnset = as.numeric(pheno$daysSinceSymptomOnset)
  pheno$ID = paste0('patient',pheno$ID)
  # sanity check
  if (!all(rownames(pheno) == colnames(exps))) {
    stop("Error: rownames of 'pheno' do not match colnames of 'exps'. Please check your data.")
  }
  
  # make a new directory for placing the result folder
  res_dir = res_dir
  if (file.exists(res_dir)){
    unlink(res_dir,recursive = TRUE)
    dir.create(res_dir,recursive = T)
  } else {
    dir.create(res_dir,recursive = T)
  }
  
  fit1 = capture.output(
    Maaslin2(input_data = exps, # capture.out is used to mutate verbose output by maaslin
             input_metadata = pheno,
             min_prevalence = 0.1, # the minimum percent of samples for which a feature is detected at minimum abundance (default 0.1)
             min_abundance = 1, # The minimum abundance for each feature (default 0)
             normalization = 'TSS', # default TSS
             transform = 'LOG', # default log 
             fixed_effects = c('daysSinceSymptomOnset'),
             #reference = 
             random_effects = "ID",
             cores = 8,
             max_significance = 0.3, # this is set to high value so that maaslin automatically plot every figure that might need in future, #FDR q value threshold
             output = res_dir))
  res = read.csv(paste(res_dir, 'all_results.tsv', sep = '/'), sep = '\t')
  return(list(exps,pheno,res))
}

# prep input matrix for maaslin. Here we just prepare the input_data, and input_metadata, and generate a directory for storing maaslin results. 
# note that we don't run maaslin2 here to allow some flexibility
prep_maaslin2_input = function(psq.tmp,feature,res_dir){
  # psq.tmp <- p.1r
  # feature <- 'Genus'
  # res_dir <- NULL
  
  exps <- taxReadCount_collapse_top10_from_phyloseq(psq.tmp,feature,sample_list = sample_names(psq.tmp))[[1]] %>%
    as.matrix %>%
    as.data.frame %>%
    setrowname
  
  exps[] <- lapply(exps, function(x) as.numeric(as.character(x)))
  if (!all(round(colSums(exps),digits = 1)==100)) {
    stop("Error: sample-wise sum of relative abundance != 100 for some samples")
  }
  # removed the unclassified features
  exps = subset(exps, !grepl('unclassified',rownames(exps)))
  
  # extract phenotype table
  pheno <- psq.tmp %>%
    sample_data %>%
    as.matrix %>%
    as.data.frame %>%
    select(c("sample.id.new","gender","race","age","vaccination","disease_severity",
             "antibiotics_antiviral_usage","treatment","BAL_Protein_Cluster",
             "daysSinceSymptomOnset","Omicron","ID","Stage","ID_Stage","Symptomonset.date")) %>%
    mutate(
      gender = as.factor(gender),
      race = as.factor(race),
      age = as.numeric(age),
      daysSinceSymptomOnset = as.numeric(daysSinceSymptomOnset),
      ID = paste0('patient', ID),
      Stage = gsub(' ','_',Stage)
    )

  # sanity check
  if (!all(rownames(pheno) == colnames(exps))) {
    stop("Error: rownames of 'pheno' do not match colnames of 'exps'. Please check your data.")
  }
  
  # make a new directory for placing the result folder
  res_dir = res_dir
  if (file.exists(res_dir)){
    unlink(res_dir,recursive = TRUE)
    dir.create(res_dir,recursive = T)
  } else {
    dir.create(res_dir,recursive = T)
  }
  return(list(exps,pheno))
}

# generate a empty directory with a given path
regenerate_dir = function(res_dir){
  if (file.exists(res_dir)){
    unlink(res_dir,recursive = TRUE)
    dir.create(res_dir,recursive = T)
  } else {
    dir.create(res_dir,recursive = T)
  }
}

# plot maaslin2 coef with error bars
plot_lfc = function(res, qcutoff,colorpalette,higher_asterisks = F){
  # test variables
  # res = res
  # qcutoff = 0.05
  
  res.sig <- res %>% filter(qval < qcutoff) %>%
    arrange(desc(coef)) %>%
    mutate(feature = as.factor(feature))
  res.sig$qval_BH_asterisks = 
    ifelse(res.sig$qval > 0.05, 'ns',
           ifelse(res.sig$qval > 0.01, "*",
                  ifelse(res.sig$qval > 0.001, '**', '***')))
  
  dodge <- position_dodge(width = 0.7)
  plt = ggplot(res.sig, aes(x = reorder(feature,coef), y = coef, fill = value, group = value)) +
    geom_bar(stat = "identity", position = dodge, width = 0.7) +
    geom_errorbar(aes(ymin = coef - stderr, ymax = coef + stderr),
                  position = dodge, width = 0.2, size = 0.4) +
    geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
    
    scale_fill_manual(values = colorpalette,name = '') +
    
    coord_flip() +
    
    labs(x = "Genus", y = "Log2 Fold Change") +
    
    theme_bw(base_size = 14) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "top",
      legend.text = element_text(size = 12, face = "bold"),
      axis.text.y = element_text(size = 12, face = "bold", color = "black"),
      axis.text.x = element_text(size = 12, color = "black"),
      axis.title = element_text(size = 14, face = "bold")
    )  
  if (higher_asterisks){
    plt = plt + geom_text(aes(y = coef + 1.8 * sign(coef), label = qval_BH_asterisks),
              position = dodge, size = 4, color = "red") 
  } else {
    plt = plt + geom_text(aes(y = coef + 1 * sign(coef), label = qval_BH_asterisks),
                          position = dodge, size = 4, color = "red") 
  }

  return(plt)
}

# summarize df by first grouping rows with "groupby" (can be a vector of multiple column names), then calculate mean, sd, se, median, q25, q75 for value column
summarize_grouped_stats <- function(df, groupby, value) {
  df %>%
    group_by(across(all_of(groupby))) %>%
    summarise(
      mean = mean(.data[[value]], na.rm = TRUE),
      sd = sd(.data[[value]], na.rm = TRUE),
      se = sd / sqrt(sum(!is.na(.data[[value]]))),
      median = median(.data[[value]], na.rm = TRUE),
      q25 = quantile(.data[[value]], 0.25, na.rm = TRUE),
      q75 = quantile(.data[[value]], 0.75, na.rm = TRUE),
      .groups = "drop"
    )
}

# save a plot to pdf
saveplot = function(plt,filename,w = 5, h = 5) {
  pdf(file.path(outdir, paste0(filename,".pdf")), width = w, height = h)
  try({print(plt)}, silent = TRUE)
  dev.off()
}
