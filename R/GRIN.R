# GRIN - Geneset Refinement using Interacting Networks

####################Description####################################################
# Filter a gene set based on distribution of ranks using Mann-Whitney U tests
# and a null distribution of ranks based on geneset size and multiplex network used.
# This command line script takes a pre-calculated interaction network, and a geneset.
# Authors: David Kainer, Kyle Sullivan, Matthew Lane, Mikaela Cashman McDevitt, Izaak Miller
###################################################################################

parse_arguments <- function() {
  suppressPackageStartupMessages(library(optparse))
  
  option_list = list(
    make_option(c("-d", "--data"),
                action="store",
                default=NULL,
                type='character',
                help="path to the .Rdata file for your combo of underlying functional networks. This file is produced by RWR_make_MHobject.R"),
    make_option(c("-g", "--geneset"), action="store", default=NULL, type='character',
                help="path to the set of seed genes. Tab-delimited file must have the following first two cols without heading: <setid> <gene>\n
                    If you wish to include weights for the genes then the third col should be numeric: <setid> <gene> <weight>\n
                    Note, the weights can be on any numeric scale (they will be normalised) but should all be > 0"),
    make_option(c("-r", "--restart"),
                action="store",
                default=0.7,
                type='numeric',
                help="set the restart parameter. Higher value means the walker will jump back to a seed node more often. default [default %default]"),
    make_option(c("-t", "--tau"),
                action="store",
                default="1.0",
                help="comma-separated list of values between that MUST add up to the number of network layers in the .Rdata file.\n 
                    One value per network layer that determines the probability that the random walker will restart in that layer.\n
                    e.g. if there are three layers (A,B,C) in your multiplex network, then --tau '0.2,1.3,1.5' will mean that layer A is 
                    less likely to be walked on after a restart than layers B or C."),
    make_option(c("-m", "--modname"),
                action="store",
                default="default",
                type='character',
                help="alias for this run. Useful for output."),
    make_option(c("-p", "--plot"),
                action="store_true",
                default=FALSE,
                help="Include this parameter if you want to output PNG plots of results. [default %default]"),
    make_option(c("-o", "--outdir"),
                action="store",
                default=NULL,
                type='character',
                help="path to the output directory"),
    make_option(c("--threads"),
                action="store",
                default=parallel::detectCores()-1,
                type='numeric',
                help="number of threads to use. default for your system is all cores - 1 [default %default]"),
    make_option(c("-s", "--simple-filenames"),
                action="store_true",
                default=FALSE,
                help="Use simple filenames."),
    make_option(c("-v", "--verbose"),
                action="store_true",
                default=FALSE,
                help="log more stuff")
  )
  
  desc <- "GRIN.R"
  opt <- parse_args(OptionParser(option_list=option_list,
                                description=desc),
                   convert_hyphens_to_underscores=TRUE)
  
  errors <- 0
  # Check whether all necessary arguments have been set by the user
  # Check opt$data.
  if(is.null(opt$data)) {
    message("ERROR:: --data is required but is not set.")
    errors <- errors+1
  } else if(!file.exists(opt$data) ) {
    message("ERROR:: --data must be an existing Rdata file.")
    errors <- errors+1
  }
  # Check opt$geneset.
  if(is.null(opt$geneset)) {
    message("ERROR:: --geneset is required but is not set.")
    errors <- errors+1
  } else if (!file.exists(opt$geneset)) {
    message("ERROR:: --geneset must be an existing TSV file.")
    errors <- errors+1
  }
  # Check outputs: opt$outdir.
  if(is.null(opt$outdir)) {
    message("ERROR:: Output directory not provided using --outdir.")
    errors <- errors+1
  }
  if(opt$plot & is.null(opt$outdir)) {
    message("ERROR:: --plot is True but --outdir was not given; --outdir is required with --plot.")
    errors <- errors+1
  }
  
  if (opt$verbose) {
    print(opt)
  }
  
  if(errors > 0) {
    quit()
  }
  
  return(opt)
}

########################################################################
# Functions
########################################################################
load_geneset = function(path, nw.mpo=NULL, opt=NULL) {
  if (is.null(path)) {
    return(NULL)
  } else {
    geneset <- read.table(path, header = F, sep="\t", colClasses=c('character'), fill=TRUE)
    if (ncol(geneset) < 2) {
      stop("Your geneset file is incorrectly formatted. See help with -h option.")
    }
    
    # Check if seed weights are included by user or not.
    if (is.numeric(geneset$V3)) {
      geneset <- dplyr::select(geneset,1:3)
      colnames(geneset) = c("setid","gene","weight")
    } else {
      # If there is a non-numeric third col, then just give all genes a weight of 1.
      geneset <- dplyr::select(geneset,1:2) %>% mutate(weight = 1)
      colnames(geneset) = c("setid","gene","weight")
    }
  }
  
  if (!is.null(opt) && opt$verbose) {
    message(sprintf('Loaded geneset (%s genes):', nrow(geneset)))
    print(head(geneset))
  }
  
  return(geneset)
}

get_or_set_tau = function(nw.mpo, opt) {
  # Ensure Tau param is appropriate (must be one value per layer, and must add up to NumLayers).
  tau <- as.numeric(unlist(strsplit(opt$tau, split=",")))
  if (sum(tau) != nw.mpo$Number_of_Layers || length(tau)!=nw.mpo$Number_of_Layers) {
    message(sprintf("WARNING:: Your comma-delimited tau parameter values do not add up to the number of network layers: %s", opt$tau))
    tau <- rep(1, nw.mpo$Number_of_Layers)
    message("WARNING:: Automatically re-setting tau = ", appendLF=FALSE)
    print(tau)
  }
  return(tau)
}

write_table = function(table, path) {
  # Create out_dir if it doesn't exist (avoid warning message if out_dir exists).
  out_dir = dirname(path)
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive=TRUE)
  }
  # Save the table.
  write.table(table, 
              path,
              sep = "\t",
              quote=F,
              col.names = T,
              row.names = F)
}


get_file_path = function(..., outdir=NULL, ext='.tsv') {
  filename = paste(..., sep='__')
  filename = paste0(filename, ext)
  if (!is.null(outdir)) {
    filename = file.path(outdir, filename)
  }
  return(filename)
}

RWR <- function(geneset, adjnorm, mpo, restart = 0.7, tau = 1, name=NULL, threads=1, verbose = NULL) {
  # Run RWR for each seed gene one at a time.
  
  # This is just the name of the combined networks.
  networks <- paste(names(mpo)[1:mpo$Number_of_Layers], collapse = "_") 
  
  # adjust verbose if output is happening inadvertently, possible set verbose = FALSE in function definition
  if (verbose) {
    message('Using geneset:')
    print(geneset)
    message("Running RWR on geneset...")
  }
  
  # Grab each gene from the geneset to use as the seed to find and run RWR using all other genes as seeds.
  doParallel::registerDoParallel(cores=threads)
  ranks <- foreach(r = 1:nrow(geneset)) %dopar% {
    leftout       <- geneset %>% dplyr::slice(r) %>% pull(gene)
    seed.genes    <- geneset %>% dplyr::filter(gene != leftout) %>% pull(gene)
    
    rwr <- Random.Walk.Restart.Multiplex(x = adjnorm, MultiplexObject = mpo,
                                         Seeds = seed.genes, r=restart,
                                         tau = tau, weights = 1 )
    #explain this pipeline
    rwr$RWRM_Results <- rwr$RWRM_Results %>%
      dplyr::mutate(rank = dplyr::min_rank(-Score)) %>%
      # for any gene that scored 0, give it worst possible rank
      dplyr::mutate(rank = if_else(Score==0, true = mpo$Number_of_Nodes_Multiplex,
                            false = rank)) %>% 
      dplyr::mutate(InGeneset = as.numeric(rwr$RWRM_Results$NodeNames %in% leftout)) %>%
      dplyr::mutate(num_in_network = nrow(geneset), left_out = leftout,
                    networks = networks, run = r, modname = name) %>%
      dplyr::relocate(left_out, .before = Score) %>%
      dplyr::relocate(modname, .before = left_out) %>%
      dplyr::select(-NodeNames)
    rwr$RWRM_Results
  }
  
  doParallel::stopImplicitCluster()
  ranks <- dplyr::bind_rows(ranks) %>% dplyr::filter(InGeneset == 1) %>%
    dplyr::arrange(rank)
  return(ranks)
}

randomGeneset <- function(genePool, n, i) {
  df <- data.frame(setid=i, gene=sample(genePool, n, replace = FALSE), weight=1)
  return(df)
}

# computed null distribution of ranks
nullDist <- function(geneset, adjnorm, mpo, restart = 0.7, tau = 1, name = NULL, threads = 1, verbose = NULL) {
  # Set random seed for deterministic output
  set.seed(42)
  
  # Do not set verbose equal to true for null distribution
  temp_verbose <- FALSE
  
  message('Generating null distribution...')

  # Only select from genes in multiplex network not contained in the input geneset
  allGenes <- as_tibble(mpo$Pool_of_Nodes)
  genePool <- unlist(allGenes %>% dplyr::filter(!value %in% geneset$V2))
  
  # Generate 100 sets of random genes for use with RWR to generate null ranks
  nullRanks <- list()
  numRand <- 100
  for (i in 1:numRand) {
    randomSet <- randomGeneset(genePool, nrow(geneset), paste0("rand",i))
    ranks.random <- RWR(randomSet, adjnorm, mpo, restart, tau, name, threads, temp_verbose)
    nullRanks[[i]] <- ranks.random
    # Generate percentages to show progress of generating null distribution
    progress <- (i/numRand)*100
    if (progress %% 5 == 0) message(sprintf('%s%% complete', progress))
  }
  
  # Only select rank in each nullRanks
  for (i in 1:length(nullRanks)) {
    nullRanks[[i]] <- nullRanks[[i]] %>% dplyr::arrange(rank) %>% dplyr::select(rank)
  }
  
  # Bind columns for all null gene ranks
  nullRanks <- dplyr::bind_cols(nullRanks, .name_repair = "unique")
  
  # Take median of ranks to generate null distribution
  nullRanks <- nullRanks %>% rowwise() %>%
    mutate(med = median(c_across(where(is.numeric)), na.rm = TRUE)) %>%
    dplyr::select(med)
  
  message('Null distribution calculated.')
  return(nullRanks)
}

# Calculate sliding window Mann-Whitney U test p-values and elbow of resulting
# curve to determine cutoff point
mannWhitneyWindow <- function(nullRanks, scores) {
  # Make sliding window of size 0.15 * number of genes
  numGenes <- nrow(scores)
  windowSize <- round(numGenes * 0.15, digits = 0)
  
  windowMatrix <- foreach(i = 1:(numGenes-windowSize), .combine = 'rbind') %do% {
    winstart  <- i
    winend    <- winstart+windowSize
    window    <- dplyr::slice(scores, winstart:winend)
    window.null <- nullRanks[winstart:winend,1]
    # Calculate Mann-Whitney U test (two-sample Wilcoxon rank sum test)
    test <- wilcox.test(window$rank, window.null$med, alternative = "less", paired = F) 
    df <- data.frame(window = i, p = test$p.value)
    colnames(df) <- c("Window", "p")
    df
  }
  return(windowMatrix)
}

elbowFilter <- function(scores, windowMatrix) {
  # Calculate elbow from sliding window and round to nearest whole number 
  elbow <- KneeArrower::findCutoff(windowMatrix$Window, windowMatrix$p, method = "first")
  elbowRound <- round(elbow$x, digits = 0)
  
  # Now filter genes into retained and removed gene sets
  retainedGenes <- dplyr::filter(scores, rank_position < elbowRound) %>%
    dplyr::mutate(set = "Retained") %>% dplyr::select(-rank_position)
  removedGenes <- dplyr::filter(scores, rank_position >= elbowRound) %>%
    dplyr::mutate(set = "Removed") %>% dplyr::select(-rank_position)
  
  filteredGenes <- list(retainedGenes, removedGenes, elbowRound)
  names(filteredGenes) <- c("Retained_Genes", "Removed_Genes", "Elbow")
  return(filteredGenes)
}

save_plot <- function(scores, windowMatrix, elbow, opt) {
  message("Generating and saving elbow plot...")
  
  # Elbow plots for ranked and randomly removed genes
  outplot <- ggplot(windowMatrix) + 
    geom_line(aes(x=Window, y=p), size=1) + 
    geom_vline(xintercept=elbow, linetype="dashed") + 
    labs(title=scores$setid[1],
         subtitle="Sliding window Mann-Whitney U Test") + 
    theme_light() + theme(axis.text.x = element_text(size=6))
    
  if (opt$simple_filenames) {
    filepath = get_file_path('GRIN-elbow-plots',opt$outdir, ext='.png')
  } else {
    filepath = get_file_path("GRIN", opt$modname, "_elbow_plot", 
                             outdir=opt$outdir, ext='.png')
  }

  png(filename=filepath, width=800, height=800)
  print(outplot)
  dev.off()
}

main <- function() {
  # Parse command line arguments
  opt <- parse_arguments()
  
  suppressPackageStartupMessages({
    library(RandomWalkRestartMH)
    library(igraph)
    library(dplyr)
    library(foreach)
    library(signal)
    library(KneeArrower)
  })
  
  # Load multiplex network and adjacency matrix for RWR
  load(opt$data)
  
  # Create output directory if it does not yet exist
  if(!is.null(opt$outdir) & !dir.exists(opt$outdir)) {
    dir.create(opt$outdir)
  }
  
  # Load gene set
  geneset.orig <- load_geneset(opt$geneset, nw.mpo, opt)
  
  # Remove any duplicate genes in the geneset
  ngenes <- nrow(geneset.orig)
  geneset <- geneset.orig %>% dplyr::distinct(gene, .keep_all=T)
  # If duplicate genes exist, write to file and remove from input for RWR
  if (nrow(geneset) < ngenes) {
    duplicates <- geneset.orig[duplicated(geneset.orig), ]
    message(sprintf('%s duplicate genes removed', nrow(duplicates)))
    print(duplicates)
    duplicates <- duplicates %>% dplyr::distinct(gene, .keep_all=T)
    filepath = get_file_path("duplicate_genes", opt$modname, outdir=opt$outdir,
                             ext='.txt')
    write_table(duplicates, filepath)
    ngenes <- nrow(geneset)
    geneset.orig <- geneset
  }
  
  # Remove any genes not in the multiplex
  geneset <- geneset.orig %>% dplyr::filter(gene %in% nw.mpo$Pool_of_Nodes)
  
  not_in_multiplex <- NULL
  # If genes not in multiplex, print genes not in multiplex and write to file
  if(nrow(geneset.orig) < ngenes) {
    message(sprintf('Only %s genes are present in the multiplex', nrow(geneset)))
    not_in_multiplex <- geneset %>%
      dplyr::slice(which(!geneset$gene %in% nw.mpo$Pool_of_Nodes))
    message(sprintf('%s genes were not found in multiplex:',nrow(not_in_multiplex)))
    message('WARNING: missing genes can affect results. Consider whether to drop genes from gene list or adjust multiplex network layers accordingly.')
    print(not_in_multiplex)
    # Write genes not found in multiplex to file
    filepath = get_file_path("genes_not_in_multiplex", opt$modname,
                             outdir=opt$outdir, ext='.txt')
    write_table(not_in_multiplex, filepath)
  } 
  else {
    message(sprintf('All %s genes in the geneset are present in the multiplex', nrow(geneset)))
  }
  
  # Throw error if geneset contains less than 4 genes
  if (nrow(geneset) < 4) {
    message('ERROR: At least 4 genes must be in the geneset for GRIN to work.')
    message('ERROR: Your geneset:')
    print(head(geneset))
    message(paste0("Number of genes in geneset: ", nrow(geneset)))
    return(0)
  }
  
  # Set tau
  tau <- get_or_set_tau(nw.mpo, opt)
  
  # Calculate null distribution
  nullDist <- nullDist(geneset, nw.adjnorm, nw.mpo, opt$restart, tau=tau,
                       opt$modname, opt$threads, opt$verbose)
  
  # Calculate RWR ranks of input gene set
  message("Filtering input gene set...")
  ranks <- RWR(geneset, nw.adjnorm, nw.mpo, opt$restart, tau, opt$modname,
               opt$threads, opt$verbose)
  
  # Break ties in rank for sliding window comparison to null distribution
  ranks <- ranks %>% dplyr::mutate(rank = rank + if_else(duplicated(rank), runif( n(),0,1 ), 0)) %>%
    dplyr::select(rank, left_out) %>% dplyr::relocate(left_out, .before = rank)

  ranks <- left_join(ranks, geneset, by = c("left_out"="gene"))
  ranks <- ranks %>% dplyr::relocate(setid, .before = left_out) %>%
    dplyr::mutate(rank_position = 0)
  # Add rank position next to each gene from input gene set
  for (i in 1:(nrow(geneset))) ranks$rank_position[i] <- i
  
  # Compute Mann-Whitney U test with sliding window
  windowMatrix <- mannWhitneyWindow(nullDist, ranks)
  
  # Compute elbow and filter genes into retained and removed gene sets
  filteredGenes <- elbowFilter(ranks, windowMatrix)
  
  # Add genes not in multiplex to removed genes
  if (!is.null(not_in_multiplex)) {
    removedGenes <- rbind(filteredGenes$Removed_Genes, not_in_multiplex, fill = TRUE)
  }
  else {
    removedGenes <- filteredGenes$Removed_Genes
  }
  
  # Write retained and removed gene sets to file
  retainedPath <- get_file_path("GRIN", opt$modname,
                                "Retained_Genes", outdir=opt$outdir, ext='.txt')
  removedPath <- get_file_path("GRIN", opt$modname,
                               "Removed_Genes", outdir=opt$outdir, ext='.txt')
  write_table(filteredGenes$Retained_Genes,retainedPath)
  write_table(removedGenes,removedPath)
  
  # If flag present, plot elbow plot and save sliding window matrix
  if(opt$plot) {
    suppressPackageStartupMessages(library(ggplot2))
    message("Making plot and saving sliding window matrix...")
    windowPath <- get_file_path("GRIN", opt$modname,
                                "Window_Matrix", outdir=opt$outdir, ext='.txt')
    write_table(windowMatrix, windowPath)
    save_plot(ranks, windowMatrix, filteredGenes$Elbow, opt)
  }

  message(paste0("COMPLETED GRIN: ",opt$modname))
  return(0)
}

status = main()
quit(save='no', status=status)
