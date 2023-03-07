#!/usr/bin/env Rscript

##################################
# create MCMC_inferCNV S4 object #
##################################
#' MCMC_inferCNV class
#'
#' @description Uses Markov Chain Monte Carlo (MCMC) and Gibbs sampling to estimate the posterior
#' probability of being in one of six Copy Number Variation states (states: 0, 0.5, 1, 1.5, 2, 3) for CNV's identified by
#' inferCNV's HMM. Posterior probabilities are found for the entire CNV cluster and each individual
#' cell line in the CNV.
#'
#' @slot bugs_model BUGS model.
#' @slot sig fitted values for cell lines, 1/standard deviation to be used for determining the distribution of each cell line
#' @slot mu Mean values to be used for determining the distribution of each cell line
#' @slot group_id ID's given to the cell clusters.
#' @slot cell_gene List containing the Cells and Genes that make up each CNV.
#' @slot cnv_probabilities Probabilities of each CNV belonging to a particular state from 0 (least likely)to 1 (most likely).
#' @slot cell_probabilities Probabilities of each cell being in a particular state, from 0 (least likely)to 1 (most likely).
#' @slot args Input arguments given by the user
#' @slot cnv_regions ID for each CNV found by the HMM
#'
#' @exportClass MCMC_inferCNV
#' @name MCMC_inferCNV-class
#' @rdname MCMC_inferCNV-class
#' @keywords classes
#'
# Requires:
# infercnv, rjags, ggplot2, parallel, futile.logger, reshape
## build off of the present S4 object inferCNV_obj to add more slots
MCMC_inferCNV <- setClass("MCMC_inferCNV", slots = c(bugs_model = "character",
                                                     sig = "numeric",
                                                     mu = "numeric",
                                                     group_id = "integer",
                                                     cell_gene = "list",
                                                     cnv_probabilities = "list",
                                                     cell_probabilities = "list",
                                                     args = "list",
                                                     cnv_regions = "factor"),
                          contains = "infercnv")




#############
# Accessors #
#############

#' Access the values for cellGene
#'
#' This function returns the list of values in cellGene
#'
#' @param obj The MCMC_inferCNV_obj S4 object.
#'
#' @return A list.
#' @rdname cellGene-method
#' @keywords internal
#' @noRd
setGeneric(name = "cellGene",
           def = function(obj) standardGeneric("cellGene"))
#' @rdname cellGene-method
#' @aliases cellGene
#' @noRd
setMethod(f = "cellGene",
          signature = "MCMC_inferCNV",
          definition=function(obj) obj@cell_gene)

#' Access the modle file
#'
#' This function returns the location of the model file 
#'
#' @param obj The MCMC_inferCNV_obj S4 object.
#'
#' @return A list.
#' @rdname modelFile-method
#' @keywords internal
#' @noRd
setGeneric(name = "modelFile",
           def = function(obj) standardGeneric("modelFile"))
#' @rdname modelFile-method
#' @aliases modelFile
#' @noRd
setMethod(f = "modelFile",
          signature = "MCMC_inferCNV",
          definition=function(obj) obj@bugs_model )

#' Access the expression data
#'
#' This function returns the expression data from object
#'
#' @param obj The MCMC_inferCNV_obj S4 object.
#'
#' @return data frame.
#' @rdname expirData-method
#' @keywords internal
#' @noRd
setGeneric(name = "expirData",
           def = function(obj) standardGeneric("expirData"))
#' @rdname expirData-method
#' @aliases expirData
#' @noRd
setMethod(f = "expirData",
          signature = "MCMC_inferCNV",
          definition=function(obj) obj@expr.data )

#' Access the arguments in the s4 object 
#'
#' Return the list of arguments passed to the inferCNVBayesNet function.
#'
#' @param obj The MCMC_inferCNV_obj S4 object.
#'
#' @return list
#' @rdname getArgs-method
#' @keywords internal
#' @noRd
setGeneric(name = "getArgs",
           def = function(obj) standardGeneric("getArgs"))
#' @rdname getArgs-method
#' @aliases getArgs
#' @noRd
setMethod(f = "getArgs",
          signature = "MCMC_inferCNV",
          definition=function(obj) obj@args)

#######################
# Object Manipulation #
#######################
#'
#' Get the cell Mean and Standard Deviation for identified cnv regions
#'
#' @param obj The MCMC_inferCNV_obj S4 object.
#' @param HMM_states HMM i3 indentified states.
#' @param infercnv_obj inferCNV object.
#'
#' @return obj The MCMC_inferCNV_obj S4 object.
#'
#' @rdname MeanSD-method
#' @keywords internal
#' @noRd
setGeneric(name="MeanSD",
           def=function(obj, HMM_states, infercnv_obj)
           { standardGeneric("MeanSD") }
)

#' @rdname MeanSD-method
#' @aliases MeanSD
#' @noRd
setMethod(f="MeanSD",
          signature="MCMC_inferCNV",
          definition=function(obj, HMM_states, infercnv_obj)
          {
              # i6 HMM method 
              if (getArgs(obj)$HMM_type == 'i6'){
                  gene_expr_by_cnv = .get_gene_expr_by_cnv(obj@.hspike)
                  cnv_mean_sd = .get_gene_expr_mean_sd_by_cnv(gene_expr_by_cnv)
                  cnv_sd <- cbind(lapply(cnv_mean_sd,function(x){x$sd}))
                  cnv_mean <- cbind(lapply(cnv_mean_sd,function(x){x$mean}))
                  ## Sort so in order of {x0,...,x1,..,x3} and get into a vector format
                  obj@mu <- unlist(cbind(cnv_mean[sort(row.names(cnv_mean)),]))
                  obj@sig <-  unlist(cbind(cnv_sd[sort(row.names(cnv_sd)),]))
                  obj@sig <- 1/(obj@sig^2)
                  if (getArgs(obj)$quietly == FALSE) {
                      print(paste("Means: ", obj@mu, collapse = ""))
                      print(paste("Sig: ", obj@sig, collapse = ""))
                  }
              } else {
                  # i3 HMM method 
                  # states <- unique(c(HMM_states))
                  # ## c(HMM_states), c(infercnv_obj@expr.data) : Converts HMM_states and expression data matrix into a vector 
                  # # MEAN
                  # means <- vapply(states, function(i) mean(c(infercnv_obj@expr.data)[which(c(HMM_states) %in% i)]), FUN.VALUE = numeric(1))
                  # names(means)<-states
                  # obj@mu <- means[sort(names(means))]
                  # # SD
                  # std <- sapply(states, function(i) sd(c(infercnv_obj@expr.data)[which(c(HMM_states) %in% i)]))
                  # names(std)<-states
                  # obj@sig <- std[sort(names(std))]
                  # obj@sig <- 1/(obj@sig^2)
                  suppressMessages(invisible(capture.output(cnv_mean_sd <- .i3HMM_get_sd_trend_by_num_cells_fit(obj))))
                  mean <- c(cnv_mean_sd$mu - cnv_mean_sd$mean_delta,
                            cnv_mean_sd$mu,
                            cnv_mean_sd$mu + cnv_mean_sd$mean_delta)
                
                  sd <- c(1/(cnv_mean_sd$sigma^2),
                          1/(cnv_mean_sd$sigma^2),
                          1/(cnv_mean_sd$sigma^2))
                  
                  obj@mu <- mean
                  obj@sig <- sd
                  
                  if (getArgs(obj)$quietly == FALSE) {
                      print(paste("Means: ", obj@mu, collapse = ""))
                      print(paste("Sig: ", obj@sig, collapse = ""))
                  }
              }
              return(obj)
          }
)

#' Add the probability threshold for the arguments in the MCMC infercnv object.
#'
#' This function adds the variable BayesMaxPNormal to the arguments slot of the the MCMC infercnv object.
#'
#' @param obj The MCMC_inferCNV_obj S4 object.
#' @param BayesMaxPNormal probability to be used as a threshold for CNV or cell removal.
#'
#' @return MCMC_inferCNV_obj S4 object.
#'
#' @rdname setBayesMaxPNormal-method
#' @keywords internal
#' @noRd
setGeneric(name = "setBayesMaxPNormal",
           def = function(obj, BayesMaxPNormal) standardGeneric("setBayesMaxPNormal"))

#' @rdname setBayesMaxPNormal-method
#' @aliases setBayesMaxPNormal
#' @noRd
setMethod(f = "setBayesMaxPNormal",
          signature = "MCMC_inferCNV",
          definition=function(obj, BayesMaxPNormal) {
              obj@args$BayesMaxPNormal <- BayesMaxPNormal
              return(obj)
              }
          )

#' Create a list that holds Genes and Cells for each separate identified CNV
#'
#' @param obj The MCMC_inferCNV_obj S4 object.
#' @param pred_cnv_genes_df Data for genes in each predicted CNV.
#' @param cell_groups_df Data for each cell in the predicted CNV's.
#'
#' @return obj The MCMC_inferCNV_obj S4 object.
#'
#' @rdname getGenesCells-method
#' @keywords internal
#' @noRd
setGeneric(name="getGenesCells",
           def=function(obj, pred_cnv_genes_df, cell_groups_df)
           { standardGeneric("getGenesCells") }
)

#' @rdname getGenesCells-method
#' @aliases getGenesCells
#' @noRd
setMethod(f="getGenesCells",
          signature="MCMC_inferCNV",
          definition=function(obj, pred_cnv_genes_df, cell_groups_df)
          {
              ## list that holds Genes and Cells for each separate identified CNV
              obj@cell_gene <- lapply(obj@cnv_regions,function(i) {
                  # subset the data to get the rows for the current CNV
                  current_cnv <- pred_cnv_genes_df[which(i == pred_cnv_genes_df$gene_region_name),]
                  # get the index for the genes that are in each cnv
                  genes <- current_cnv$gene
                  # pred_cnv_genes_df[which(pred_cnv_genes_df$gene_region_name %in% i),]$gene
                  gene_idx <- which(row.names(obj@expr.data) %in% genes)

                  # get the index for the cells that are in each cnv
                  sub_cells <- unique(current_cnv$cell_group_name)
                  cells_idx <- which(colnames(obj@expr.data) %in% cell_groups_df[which(cell_groups_df$cell_group_name %in% sub_cells),]$cell)

                  # get unique current CNV state
                  state <- unique(current_cnv$state)
                  return(list("cnv_regions" = i, "Genes" = gene_idx, "Cells" = cells_idx, "State" = state))
              })
              return(obj)
          }
)




#' Initialize the MCMC_inferCNV_obj object
#'
#' @param obj The MCMC_inferCNV_obj S4 object.
#' @param args_parsed The arguments given to the function.
#' @param infercnv_obj InferCNV object.
#'
#' @return obj The MCMC_inferCNV_obj S4 object.
#'
#' @rdname initializeObject-method
#' @keywords internal
#' @noRd
setGeneric(name="initializeObject",
           def=function(obj, args_parsed, infercnv_obj)
           { standardGeneric("initializeObject") }
)

#' @rdname initializeObject-method
#' @aliases initializeObject
#' @noRd
setMethod(f="initializeObject",
          signature="MCMC_inferCNV",
          definition=function(obj, args_parsed, infercnv_obj)
    {
        futile.logger::flog.info(paste("Initializing new MCM InferCNV Object."))
        files <- list.files(args_parsed$file_dir, full.names = TRUE)
  
        # Validate the inferCNV Object
        validate_infercnv_obj(infercnv_obj)
  
        ## create the S4 object
        obj <- MCMC_inferCNV(infercnv_obj)
        ## add the command line arguments
        obj@args <- args_parsed
  
        ## Load the files for cnv predictions
        cell_groups_PATH <- files[grep(files, pattern = paste0("19_HMM_pred", args_parsed$resume_file_token, ".cell_groupings"))]
        pred_cnv_genes_PATH <- files[grep(files, pattern = paste0("19_HMM_pred", args_parsed$resume_file_token, ".pred_cnv_genes.dat"))]
        cell_groups_df <- read.table(cell_groups_PATH, header = T, check.names = FALSE, sep="\t")
        pred_cnv_genes_df <- read.table(pred_cnv_genes_PATH, header = T, check.names = FALSE, sep="\t", stringsAsFactors = TRUE)
  
        # cnv region id's
        obj@cnv_regions <- unique(pred_cnv_genes_df$gene_region_name)
        futile.logger::flog.info(paste("Total CNV's: ", length(obj@cnv_regions)))
  
        ## Load Mixture Model File
        futile.logger::flog.info(paste("Loading BUGS Model."))
        obj@bugs_model <- getArgs(obj)$model_file
  
        ## list that holds Genes and Cells for each separate identified CNV
        obj <- getGenesCells(obj, pred_cnv_genes_df, cell_groups_df)
  
        # Create numerical ids for each subgroup of cells
        ## group name ids
        cell_group_id <- unique(pred_cnv_genes_df$cell_group_name)
  
        ## set numerical id's for cell groups and set values in a vector for cell positions in the matrix
        group_id <- rep(NA, max(unlist(obj@observation_grouped_cell_indices)))
        lapply(seq_along(cell_group_id), function(i) {
            ## cells in the cluster group
            cells <- cell_groups_df[cell_groups_df$cell_group_name %in% cell_group_id[i],]$cell
            ## set the numerical id in the vector
            group_id[which(colnames(obj@expr.data) %in% cells)] <<- i
        })
        obj@group_id <- group_id
  
        return(obj)
    }
)



#' Set the probabilities for each CNV belonging to each state as well as probability of each cell belonging to a states
#'
#' @param obj The MCMC_inferCNV_obj S4 object.
#' @param mcmc Sampling data.
#'
#' @return obj The MCMC_inferCNV_obj S4 object.
#'
#' @rdname getProbabilities-method
#' @keywords internal
#' @noRd
setGeneric(name="getProbabilities",
           def=function(obj, mcmc)
           { standardGeneric("getProbabilities") }
)

#' @rdname getProbabilities-method
#' @aliases getProbabilities
#' @noRd
setMethod(f="getProbabilities",
          signature="MCMC_inferCNV",
          definition=function(obj, mcmc)
          {
              ## List holding state probabilities for each CNV
              cnv_probabilities <- list()
              ## List for combining the chains in each simulation
              combined_mcmc <- list()
              ## list holding the frequency of epsilon values for each cell line
              ##  for each cnv region and subgroup
              cell_probabilities <- list()

              combinedMCMC <-
                  for(j in seq_along(mcmc)){
                      # combine the chains
                      combined_mcmc[[j]] <- do.call(rbind, mcmc[[j]])
                      # run function to get probabilities
                      ## Thetas
                      cnv_probabilities[[j]] <- cnv_prob(combined_mcmc[[j]])
                      ## Epsilons
                      cell_probabilities[[j]] <- cell_prob(combined_mcmc[[j]], obj)
                  }
              obj@cnv_probabilities <- cnv_probabilities
              obj@cell_probabilities <- cell_probabilities
              return(obj)
          }
)

#' Run simulations in Parallel
#'
#' @param obj The MCMC_inferCNV_obj S4 object.
#'
#' @return Sampling data.
#'
#' @rdname withParallel-method
#' @keywords internal
#' @noRd
setGeneric(name="withParallel",
           def=function(obj)
           { standardGeneric("withParallel") }
)

#' @rdname withParallel-method
#' @aliases withParallel
#' @noRd
setMethod(f="withParallel",
          signature="MCMC_inferCNV",
          definition=function(obj)
          {
              par_func <- function(i){
                  if (getArgs(obj)$quietly == FALSE) {
                      futile.logger::flog.info(paste("Sampling Number: ", i))
                  }
                  if(!(length(obj@cell_gene[[i]]$Cells) == 0)){
                      tumor_grouping <- obj@group_id[ obj@cell_gene[[i]]$Cells ] # subset the tumor ids for the cells wanted
                      gene_exp <- obj@expr.data[obj@cell_gene[[i]]$Genes, obj@cell_gene[[i]]$Cells]
                      return(run_gibb_sampling(gene_exp, obj))
                  } else {
                      return(list(NULL))
                  }
              }
              mc.cores = ifelse(.Platform$OS.type == 'unix', as.integer(getArgs(obj)$CORES), 1) # if windows, can only use 1 here
              futile.logger::flog.info(paste("Running Sampling Using Parallel with ", getArgs(obj)$CORES, "Cores"))
              mcmc <- parallel::mclapply(seq_along(obj@cell_gene),
                                         FUN = par_func,
                                         mc.cores = mc.cores)
              return(mcmc)
          }
)

#' Run simulations in Non-Parallel mode
#'
#' @param obj The MCMC_inferCNV_obj S4 object.
#'
#' @return Sampling data
#'
#' @rdname nonParallel-method
#' @keywords internal
#' @noRd
setGeneric(name="nonParallel",
           def=function(obj)
           { standardGeneric("nonParallel") }
)

#' @rdname nonParallel-method
#' @aliases nonParallel
#' @noRd
setMethod(f="nonParallel",
          signature="MCMC_inferCNV",
          definition=function(obj)
          {
              futile.logger::flog.info(paste("Running Gibbs sampling in Non-Parallel Mode."))
              # Iterate over the CNV's and run the Gibbs sampling.
              mcmc <- lapply(seq_along(obj@cell_gene), function(i){
                  if (getArgs(obj)$quietly == FALSE) {
                      futile.logger::flog.info(paste("Sample Number: ", i))
                  }
                  if(!(length(obj@cell_gene[[i]]$Cells) == 0)){
                      tumor_grouping <- obj@group_id[ obj@cell_gene[[i]]$Cells ] # subset the tumor ids for the cells wanted
                      gene_exp <- obj@expr.data[obj@cell_gene[[i]]$Genes, obj@cell_gene[[i]]$Cells]
                      return(run_gibb_sampling(gene_exp, obj))
                  } else {
                      return(list(NULL))
                  }
              })
              return(mcmc)
          }
)


#' Given the CNV associated probability of belonging to each possible state, 
#' reassign the state assignments made by the HMM to the state that has the highest probability. 
#'
#' @param obj The MCMC_inferCNV_obj S4 object.
#' @param HMM_states InferCNV object with HMM states in expression data.
#'
#' @return obj The MCMC_inferCNV_obj S4 object.
#'
#' @rdname reassignCNV-method
#' @keywords internal
#' @noRd
setGeneric(name="reassignCNV",
           def=function(obj, HMM_states)
           { standardGeneric("reassignCNV") }
)

#' @rdname reassignCNV-method
#' @aliases reassignCNV
#' @noRd
setMethod(f="reassignCNV",
          signature="MCMC_inferCNV",
          definition=function(obj, HMM_states)
          {
              ## reassignCNV(obj, HMM_states)
              ## reassign the state assignments made by the HMM to the state that has the highest probability 
              ##    1. identify the CNVs that have a higher probability of being in a state different than what was assigned to the CNV by the HMM. 
              ##    2. if normal state is the highest probability, ignore because this is handled by the BayesMaxPNormal threshold. 
              ##    3. Reassign the states to the new states in the HMM identified state matrix. 

              futile.logger::flog.info("Reassigning CNVs based on state probabilities.")
              # Assign state that represents normal based on the HMM method 
              normalID <- ifelse(getArgs(obj)$HMM_type == 'i6', 3, 2)

              # 1.
              # get probabilities for each cnv
              cnv_means <- sapply(obj@cnv_probabilities,function(i) colMeans(i))
              # get HMM state asignments per CNV 
              hmm_assigned_states <- vapply(cellGene(obj), function(i){i$State}, FUN.VALUE = numeric(1) )
              # get the highest probabillity state for each cnv 
              highest_prob <- apply(cnv_means, 2, function(i) which( i == max(i)))

              # 2. 
              # Handle cnvs that have normal state as its highest probability
              #     Keep original state assignment by HMM 
              # normIDX <- which(highest_prob == normalID)
              # highest_prob[normIDX] <- hmm_assigned_states[normIDX]
              
              # Change the state assignment in objects gene_cell information 
              lapply(seq_len(length(cellGene(obj))), function(i) obj@cell_gene[[i]]$State <<- highest_prob[i])

              # 3. 
              lapply(cellGene(obj), function(i){
                  HMM_states[i$Genes , i$Cells ] <<- i$State
              })
              reassignIDX <- which(hmm_assigned_states != highest_prob)

              ## If any CNVs are reassigned
              ##    Send a messge of which CNVs are being reassigned from what state -> to what new state 
              if (length(reassignIDX) > 0){
                  cnvMessage <- sapply(reassignIDX, function(i) {
                      temp <- obj@cell_gene[[i]]
                      paste(temp$cnv_regions,":", hmm_assigned_states[i], " (P=",cnv_means[hmm_assigned_states[i],i],") -> ", highest_prob[i], "(P=",cnv_means[highest_prob[i],i],")")
                  })
                  futile.logger::flog.info(paste("Changing the following CNV's states assigned by the HMM to the following based on the CNV's state probabilities.\n", paste(cnvMessage, sep = "",collapse = "\n")))
              }

              return(list(obj, HMM_states))
          }
)


#' Run simulations and remove CNV's that have a probability of being normal above a set thresholld.
#' This removes possible false posotives identified by the HMM.
#'
#' @param obj The MCMC_inferCNV_obj S4 object.
#' @param HMM_states InferCNV object with HMM states in expression data.
#'
#' @return obj The MCMC_inferCNV_obj S4 object.
#'
#' @rdname removeCNV-method
#' @keywords internal
#' @noRd
setGeneric(name="removeCNV",
           def=function(obj, HMM_states)
           { standardGeneric("removeCNV") }
)

#' @rdname removeCNV-method
#' @aliases removeCNV
#' @noRd
setMethod(f="removeCNV",
          signature="MCMC_inferCNV",
          definition=function(obj, HMM_states)
          {
              ## removeCNV(obj, HMM_states)
              ## If there are CNVs that have probabilities of being normal state greater than the applied threshold, 
              ##    1. find which CNVs to remove, 
              ##    2. reassign there state to the normal state in the state matrix 
              ##    3. remove the CNVs form; list of CNVs (cell_gene), CNV probabilities, cell probabilities 
              ##    4. Write the state probabilities for each CNV to a table.
              
              # Assign index and state that represents normal based on the HMM method 
              normalID <- ifelse(getArgs(obj)$HMM_type == 'i6', 3, 2)
              # Mean values of the probability distribution of the CNV states p(CNV == {states 1:6})
              cnv_means <- sapply(obj@cnv_probabilities,function(i) colMeans(i))
              
              futile.logger::flog.info(paste("Attempting to removing CNV(s) with a probability of being normal above ", getArgs(obj)$BayesMaxPNormal))
              futile.logger::flog.info(paste("Removing ",length(which(cnv_means[normalID,] > getArgs(obj)$BayesMaxPNormal)), " CNV(s) identified by the HMM."))
              
              # If no CNV's need to be removed, stop running function and return the object and HMM_states 
              if (length(which(cnv_means[normalID,] > getArgs(obj)$BayesMaxPNormal)) == 0){ return(list(obj, HMM_states)) }
              
              # check if any CNVs with probability of being normal greater than threshold
              if (any(cnv_means[normalID,] > getArgs(obj)$BayesMaxPNormal)){

                  # 1.
                  remove_cnv <- which(cnv_means[normalID,] > getArgs(obj)$BayesMaxPNormal)
                  
                  if (getArgs(obj)$quietly == FALSE) { print("CNV's being removed have the following posterior probabilities of being a normal state: ") }
                  
                  # 2. 
                  lapply(remove_cnv, function(i) {
                      if (getArgs(obj)$quietly == FALSE) {
                          print( paste(cellGene(obj)[[i]]$cnv_regions, ", Genes: ", length(cellGene(obj)[[i]]$Genes), " Cells: ", length(cellGene(obj)[[i]]$Cells)) )
                          # print(paste(paste( "Probabilities: "), cnv_means[,i]))
                      }
                      ## Change the states to normal states
                      HMM_states[cellGene(obj)[[i]]$Genes , cellGene(obj)[[i]]$Cells ] <<- normalID
                  })

                  # 3. 
                  ## Remove the CNV's from the following matrices
                  obj@cell_gene <- obj@cell_gene[-remove_cnv]
                  obj@cell_probabilities <- obj@cell_probabilities[-remove_cnv]
                  obj@cnv_probabilities <- obj@cnv_probabilities[-remove_cnv]
                  cnv_means <- cnv_means[,-remove_cnv]
                  # obj@mcmc <- obj@mcmc[-remove_cnv]
                  # obj@combined_mcmc <- obj@combined_mcmc[-remove_cnv]
                  futile.logger::flog.info(paste("Total CNV's after removing: ", length(obj@cell_gene)))
              }
              
              # 4.
              # Write the state probabilities for each CNV to a table.
              ## check if output directory exists, if not create it 
              if(getArgs(obj)$out_dir != "." & !file.exists(getArgs(obj)$out_dir)){
                  # create the output directory
                  dir.create(file.path(getArgs(obj)$out_dir))
                  futile.logger::flog.info(paste("Creating the following Directory: ", getArgs(obj)$out_dir))
              }
              ## set column names to the CNV ID
              cnv_regions <- sapply(obj@cell_gene, function(i) { as.character(i$cnv_regions) })
              colnames(cnv_means) <- cnv_regions
              ## set row names to the states 1:6 or 1:3
              temp <- ifelse(getArgs(obj)$HMM_type == 'i6', 6, 3)
              row.names(cnv_means) <- c(sprintf("State:%s",seq_len(temp)))
              write.table(cnv_means,file = file.path(getArgs(obj)$out_dir, "CNV_State_Probabilities.dat"), col.names = TRUE, row.names=TRUE, quote=FALSE, sep="\t")
              return(list(obj, HMM_states))
          }
)

#' Run simulations and remove cells from cnv's that are predicted to be normal
#'
#' @param obj The MCMC_inferCNV_obj S4 object.
#' @param HMM_states InferCNV object with HMM states in expression data.
#'
#' @return obj The MCMC_inferCNV_obj S4 object.
#'
#' @rdname removeCells-method
#' @keywords internal
#' @noRd
setGeneric(name="removeCells",
           def=function(obj, HMM_states)
           { standardGeneric("removeCells") }
)

#' @rdname removeCells-method
#' @aliases removeCells
#' @noRd
setMethod(f="removeCells",
          signature="MCMC_inferCNV",
          definition=function(obj, HMM_states)
          {
              if (getArgs(obj)$HMM_type == 'i6'){
                  if (any(do.call(cbind, obj@cell_probabilities)[3,] > getArgs(obj)$BayesMaxPNormal)){
                      lapply(seq_along(obj@cell_probabilities), function(i) {
                          idx <- which(obj@cell_probabilities[[i]][3,] > getArgs(obj)$BayesMaxPNormal)
                          if(length(idx) > 0){
                              ## change the states to normal states
                              HMM_states[ obj@cell_gene[[i]]$Genes , obj@cell_gene[[i]]$Cells[ idx ] ] <<- 3
                              ## remove these cells from the cnv
                              obj@cell_gene[[i]]$Cells <<- obj@cell_gene[[i]]$Cells[- idx]
                          }
                      })
                      # recursively run again
                      obj <- runMCMC(obj)
                  }
              } else {
                  if (any(do.call(cbind, obj@cell_probabilities)[2,] > getArgs(obj)$BayesMaxPNormal)){
                      lapply(seq_along(obj@cell_probabilities), function(i) {
                          idx <- which(obj@cell_probabilities[[i]][2,] > getArgs(obj)$BayesMaxPNormal)
                          if(length(idx) > 0){
                              ## change the states to normal states
                              HMM_states[ obj@cell_gene[[i]]$Genes , obj@cell_gene[[i]]$Cells[ idx ] ] <<- 2
                              ## remove these cells from the cnv
                              obj@cell_gene[[i]]$Cells <<- obj@cell_gene[[i]]$Cells[- idx]
                          }
                      })
                      # recursively run again
                      obj <- runMCMC(obj)
                  }
              }
              return(obj)
          }
)

#' Run simulations using rjags.
#'
#' Run MCMC simulations using rjags. Also returns a plot the probability of each CNV being
#' normal before running any kind of post MCMC modification.
#'
#' @param obj The MCMC_inferCNV_obj S4 object.
#' 
#' @param diagnostics Option to create the diagnostic plots. 
#'
#' @return obj The MCMC_inferCNV_obj S4 object.
#'
#' @rdname runMCMC-method
#' @keywords internal
#' @noRd
setGeneric(name="runMCMC",
           def=function(obj, diagnostics)
           { standardGeneric("runMCMC") }
)

#' @rdname runMCMC-method
#' @aliases runMCMC
#' @noRd
setMethod(f="runMCMC",
          signature="MCMC_inferCNV",
          definition=function(obj,
                              diagnostics = FALSE)
          {
              # Run MCMC
              if(getArgs(obj)$CORES == 1){
                  mcmc <- nonParallel(obj)
              } else {
                  mcmc <- withParallel(obj)
              }
              
              # Create Diagnostic Plots.
              if (diagnostics == TRUE){
                  mcmcDiagnosticPlots(obj, mcmc)
              }

              # Get the probability of of each cell line and complete CNV belonging to a specific state
              futile.logger::flog.info(paste("Obtaining probabilities post-sampling"))
              obj <- getProbabilities(obj,mcmc)
              return(obj)
          }
)


##########################
# Plotting functions #
##########################
#' Get the probability of each cnv being a normal state and plot these probabilities.
#'
#' @param obj The MCMC_inferCNV_obj S4 object.
#' @param PNormal Option to add specific title to plot.
#' @param title Title to be used for the plot.
#' @param output_filename The file to write to. 
#'
#' @return obj The MCMC_inferCNV_obj S4 object.
#'
#' @rdname postProbNormal-method
#' @keywords internal
#' @noRd
setGeneric(name="postProbNormal",
           def=function(obj, PNormal, title, output_filename, useRaster)
           { standardGeneric("postProbNormal") }
)

#' @rdname postProbNormal-method
#' @aliases postProbNormal
#' @noRd
setMethod(f="postProbNormal",
          signature="MCMC_inferCNV",
          definition=function(obj, PNormal, title, output_filename, useRaster)
          {
              if (getArgs(obj)$plotingProbs == TRUE){
                  # get probability of the cnv's belonging to each state
                  cnv_means <- sapply(obj@cnv_probabilities,function(i) colMeans(i))
                  # Adjust the probabilities so greater probability corresponds to less likely to be normal
                  if ( getArgs(obj)$HMM_type == 'i6'){
                      normal_prob <- 1 - cnv_means[3,]
                  } else {
                      normal_prob <- 1 - cnv_means[2,]
                  }
                  obj@expr.data[,] <- 0
                  lapply(seq_along(normal_prob), function(i) {
                      ## change the states to normal states
                      obj@expr.data[obj@cell_gene[[i]]$Genes , obj@cell_gene[[i]]$Cells ] <<- normal_prob[i]
                  })
                  infercnv::plot_cnv(infercnv_obj          = obj,
                                     out_dir               = getArgs(obj)$out_dir,
                                     k_obs_groups          = getArgs(obj)$k_obs_groups,
                                     cluster_by_groups     = getArgs(obj)$cluster_by_groups,
                                     title                 = title,
                                     output_filename       = output_filename,
                                     write_expr_matrix     = FALSE,
                                     x.center              = 0,
                                     x.range               = c(0,1),
                                     useRaster             = useRaster
                  )
              }
          }
)

#' Plots the probability for each cnv belonging to a specific state and the probability of
#' each cell line belonging to a specific states.
#'
#' @param obj The MCMC_inferCNV_obj S4 object.
#'
#' @return obj The MCMC_inferCNV_obj S4 object.
#'
#' @rdname plotProbabilities-method
#' @keywords internal
#' @noRd
setGeneric(name="plotProbabilities",
           def=function(obj)
           { standardGeneric("plotProbabilities") }
)

#' @rdname plotProbabilities-method
#' @aliases plotProbabilities
#' @noRd
setMethod(f="plotProbabilities",
          signature="MCMC_inferCNV",
          definition=function(obj)
          {
              if (getArgs(obj)$plotingProbs == TRUE){
                  futile.logger::flog.info(paste("Creating Plots for CNV and cell Probabilities."))
                  # Plotting
                  ## plots the probability of each cell line being a particular state
                  ## plots the probability of a cnv being a particular state

                  ## add threshold to the plot title if given
                  if (!is.null(getArgs(obj)$BayesMaxPNormal)) {
                      file_CELLplot <- sprintf("cellProbs.%s.pdf",getArgs(obj)$BayesMaxPNormal)
                  } else{
                      file_CELLplot <- "cellProbs.pdf"
                  }
                  pdf(file = file.path(file.path(getArgs(obj)$out_dir),file_CELLplot), onefile = TRUE)
                  lapply(seq_along(obj@cell_probabilities), function(i){
                      print(plot_cell_prob(as.data.frame(obj@cell_probabilities[[i]]), as.character(obj@cell_gene[[i]]$cnv_regions), getArgs(obj)$HMM_type))
                  })
                  dev.off()

                  ## Plot the probability of each state for a CNV
                  ## add threshold to the plot title if given
                  if (!is.null(getArgs(obj)$BayesMaxPNormal)) {
                      file_CNVplot <- sprintf("cnvProbs.%s.pdf",getArgs(obj)$BayesMaxPNormal)
                  } else{
                      file_CNVplot <- "cnvProbs.pdf"
                  }
                  pdf(file = file.path(file.path(getArgs(obj)$out_dir), file_CNVplot), onefile = TRUE)
                  lapply(seq_along(obj@cell_probabilities), function(i){
                      print(plot_cnv_prob(obj@cnv_probabilities[[i]], as.character(obj@cell_gene[[i]]$cnv_regions), getArgs(obj)$HMM_type))
                  })
                  dev.off()
              }
          }
)

#' Create Diagnostic Plots And Summaries.
#'
#' Create Diagnostic Plots And Summaries in order to determine if convergence has occured.
#'
#' @param obj The MCMC_inferCNV_obj S4 object.
#' @param mcmc Sampling data.
#'
#' @return obj The MCMC_inferCNV_obj S4 object.
#'
#' @rdname mcmcDiagnosticPlots-method
#' @keywords internal
#' @noRd
setGeneric(name="mcmcDiagnosticPlots",
           def=function(obj, mcmc)
           { standardGeneric("mcmcDiagnosticPlots") }
)

#' @rdname mcmcDiagnosticPlots-method
#' @aliases mcmcDiagnosticPlots
#' @noRd
setMethod(f="mcmcDiagnosticPlots",
          signature="MCMC_inferCNV",
          definition=function(obj, mcmc)
          {
              futile.logger::flog.info(paste("Creating Diagnostic Plots."))
              ###########################
              # trace and denisty plots
              ###########################
              #--------------------------------------
              # trace and denisty plots for each cnv
              #--------------------------------------
              ## get the theta values
              if (getArgs(obj)$quietly == FALSE) { futile.logger::flog.info(paste("Plotting CNV Trace and Density Plots.")) }
              cnvProb <- function(combined_samples) {
                  thetas <- combined_samples[,grepl('theta', colnames(combined_samples))]
              }
              cnvMCMCList <- lapply(seq_along(mcmc), function(i){
                  lapply(mcmc[[i]], cnvProb)
              })
              # trace and denisty plots
              pdf(file = file.path(file.path(getArgs(obj)$out_dir),"CNVDiagnosticPlots.pdf"), onefile = TRUE)
              lapply(seq_along(cnvMCMCList), function(i){
                  plot(coda::mcmc.list(cnvMCMCList[[i]]))
              })
              dev.off()

              #---------------------------------------
              # trace and denisty plots for each cell
              #---------------------------------------
              ## get the theta values
              if (getArgs(obj)$quietly == FALSE) { futile.logger::flog.info(paste("Plotting Cell Trace and Density Plots.")) }
              cellProb <- function(samples) {
                  epsilons <- samples[,grepl('epsilon', colnames(samples))]
              }

              cellMCMCList <- lapply(seq_along(mcmc), function(i){
                  lapply(mcmc[[i]], cellProb)
              })
              # trace and denisty plots
              pdf(file = file.path(file.path(getArgs(obj)$out_dir),"CellDiagnosticPlots.pdf"), onefile = TRUE)
              lapply(seq_along(cellMCMCList), function(i){
                  plot(coda::mcmc.list(cellMCMCList[[i]]))
              })
              dev.off()


              ###########################
              # Auto Correlation Plots
              ###########################
              #---------------------------------------
              # Auto Correlation for each CNV
              #---------------------------------------
              if (getArgs(obj)$quietly == FALSE) { futile.logger::flog.info(paste("Plotting CNV Autocorrelation Plots.")) }
              pdf(file = file.path(file.path(getArgs(obj)$out_dir),"CNVautocorrelationPlots.pdf"), onefile = TRUE)
              lapply(seq_along(cnvMCMCList), function(i){
                  coda::autocorr.plot(coda::mcmc.list(cnvMCMCList[[i]]))
              })
              dev.off()

              ###########################
              # Gelman Plots
              ###########################
              #---------------------------------------
              # Gelman for each CNV
              #---------------------------------------
              if (getArgs(obj)$quietly == FALSE) { futile.logger::flog.info(paste("Plotting CNV Gelman Plots.")) }
              pdf(file = file.path(file.path(getArgs(obj)$out_dir),"CNVGelmanPlots.pdf"), onefile = TRUE)
              lapply(seq_along(cellMCMCList), function(i){
                  coda::gelman.plot(coda::mcmc.list(cnvMCMCList[[i]]))
              })
              dev.off()

              ###########################
              # Summary Tables
              ###########################
              if (getArgs(obj)$quietly == FALSE) { futile.logger::flog.info(paste("Creating CNV Statistical Summary Tables.")) }
              # Function to initialize the summary tables
              theta_table <- function(x,y,w){
                  mu<- unlist(summary(x[[1]][,w])[[1]][,1])
                  stdev<- unlist(summary(x[[1]][,w])[[1]][,2])
                  q2.5<- unlist(summary(x[[1]][,w])[[2]][,1])
                  q50<- unlist(summary(x[[1]][,w])[[2]][,3])
                  q97.5<- unlist(summary(x[[1]][,w])[[2]][,5])
                  gewek = unlist(geweke.diag(x[[1]][,w], frac1=0.1, frac2=0.5))[seq_along(w)]
                  df = data.frame(mu,stdev,q2.5,q50,q97.5,gewek)
                  colnames(df) <- c('Mean','St.Dev','2.5%','50%','97.5%', "Geweke")
                  rownames(df) <- c(w)
                  #return(knitr::kable(df, caption = y))
                  return(df)
              }
              # Function to get the theta (state CNV probabilities) values
              getThetas <- function(df){ df[,grepl('theta', colnames(df))] }
              # List of statistical summary tables
              summary_table <- lapply(seq_along(mcmc), function(i) {
                  title <- sprintf("CNV %s Summary Table", obj@cell_gene[[i]]$cnv_regions)
                  thetas <- lapply(mcmc[[i]], function(x) getThetas(x))
                  w = row.names(summary(as.mcmc(thetas))[[1]])
                  return(theta_table(coda::as.mcmc(thetas), title, w))
              })
              # Theme for the grob tables
              theme.1 <- gridExtra::ttheme_default(core = list(fg_params = list(parse=TRUE, cex = 0.5)),
                                                   colhead = list(fg_params=list(parse=TRUE, cex = 0.5)),
                                                   rowhead = list(fg_params=list(parse=TRUE, cex = 0.5)))
              # List of tables, table for each CNV
              plot_list <- lapply(seq_along(summary_table), function(i) {
                  ## Create table grob object
                  table <- gridExtra::tableGrob(summary_table[[i]],rows = c("State 1","State 2","State 3","State 4","State 5","State 6"), theme = theme.1)
                  ## Create the title for the table as a seperate grob object
                  title <- sprintf("%s CNV Summary Table", obj@cell_gene[[i]]$cnv_regions)
                  title <- gridExtra::tableGrob(summary_table[[i]][1,1],rows=NULL, cols=c(title))
                  ## Combine the summary table grob and the title grob
                  tab <- gridExtra::gtable_combine(title[1,], table, along=2)
                  # Adjust the position of the title
                  tab$layout[1, c("l", "r")] <- c(7, 2)
                  tab$layout[2, c("l", "r")] <- c(7, 2)
                  return(tab)
              })
              # Combine all the tablles together as one column
              test <- gridExtra::gtable_combine(plot_list, along = 2)
              # Save the tables to a PDF document
              pdf(file = file.path(file.path(getArgs(obj)$out_dir),"CNVSummaryTablels.pdf") , paper = "a4", onefile = TRUE, height = 0, width = 0)
              print(gridExtra::marrangeGrob(grobs = test, nrow = 5, ncol = 1))
              dev.off()
          }
)



##########################
# Command line arguments #
##########################
# pargs <- optparse::OptionParser()
pargs <- argparse::ArgumentParser()
pargs$add_argument(c("-f", "--infercnv_dir"),
                   type="character",
                   action='store_true',
                   dest="file_dir",
                   metavar="File_Directory",
                   help=paste("Path to files created by inferCNV.",
                              "[Default %default][REQUIRED]"))
pargs$add_argument(c("-m", "--model"),
                   type="character",
                   action='store_true',
                   dest="model_file",
                   metavar="Model_File_Path",
                   help=paste("Path to the BUGS Model file.",
                              "[Default %default][REQUIRED]"))
pargs$add_argument(c("-p","--parallel"),
                   type="character",
                   action='store_true',
                   dest="CORES",
                   default = NULL,
                   metavar="Number_of_Cores",
                   help=paste("Option to run parallel by specifying the number of cores to be used.",
                              "[Default %default]"))
pargs$add_argument(c("-o","--out_dir"),
                   type="character",
                   action='store_true',
                   dest="out_dir",
                   default = NULL,
                   metavar="Output_Directory",
                   help=paste("Option to set the output directory to save the outputs.",
                              "[Default %default]"))
pargs$add_argument(c("-M","--method"),
                   type="character",
                   action='store_true',
                   dest="postMcmcMethod",
                   default = NULL,
                   metavar="Posterior_MCMC_Method",
                   help=paste("What actions to take after finishing the MCMC.",
                              "[Default %default]"))
pargs$add_argument(c("-x","--plot"),
                   type="logical",
                   action='store_true',
                   dest="plotingProbs",
                   default = TRUE,
                   metavar="Plot_Probabilities",
                   help=paste("Plot the posterior probabilites for each CNV and each cell line in each cnv.",
                              "[Default %default]"))

# Function to run the mixture model for given expression data
# Runs on each cnv seperately. Then seperate by tumor subgroup.
#       input:
#               gene_exp            : Gene expression data
#               MCMC_inferCNV_obj   : MCMC_inferCNV object
#       return:
#               samples     : Results of the sampleing process
#
run_gibb_sampling <- function(gene_exp,
                              MCMC_inferCNV_obj
                              ){
    if (is.null(ncol(gene_exp))){
        gene_exp <- data.frame(gene_exp)
    }
    C = ncol(gene_exp)
    G = nrow(gene_exp)
    if (getArgs(MCMC_inferCNV_obj)$quietly == FALSE) {
        futile.logger::flog.info(paste("Cells: ",C))
        futile.logger::flog.info(paste("Genes: ",G))
    }
    # quiet=FALSE
    # make Data list for model
    data <- list(
        'C' = C,                                 # number of cell lines
        'G' = G,                                 # number of genes
        'gexp' = gene_exp,                       # expression data
        'sig' = MCMC_inferCNV_obj@sig,           # fitted values for cell lines, 1/standard deviation to be used for determining the distribution of each cell line
        'mu' = MCMC_inferCNV_obj@mu              # Mean values to be used for determining the distribution of each cell line
    )
    # set initial values for each cell line begining states
    if (getArgs(MCMC_inferCNV_obj)$HMM_type == "i6"){
        # i6 method 
        inits <- list(
            list(epsilon = rep(1, C)),
            list(epsilon = rep(2, C)),
            list(epsilon = rep(3, C)),
            list(epsilon = rep(4, C)),
            list(epsilon = rep(5, C)),
            list(epsilon = rep(6, C))
        )
    } else {
        # i3 method 
        inits <- list(
            list(epsilon = rep(1, C)),
            list(epsilon = rep(2, C)),
            list(epsilon = rep(3, C))
        )   
    }
    # Create the model for rjags
    model <- rjags::jags.model(modelFile(MCMC_inferCNV_obj),
                               data     = data,
                               inits    = inits, # (Initialization) optional specification of initial values in the form of a list or a function
                               n.chains = ifelse(getArgs(MCMC_inferCNV_obj)$HMM_type == "i3", 3, 6),  # the number of parallel chains for the model
                               n.adapt  = 500, # the number of iterations for adaptation (burn in)
                               quiet    = getArgs(MCMC_inferCNV_obj)$quietly)
    stats::update(model, 200, progress.bar=ifelse(getArgs(MCMC_inferCNV_obj)$quietly,"none","text"))
    # run the rjags model
    ## set the parameters to return from sampling
    parameters <- c('theta', 'epsilon')
    samples <- rjags::coda.samples(model, parameters, n.iter=1000, progress.bar=ifelse(getArgs(MCMC_inferCNV_obj)$quietly,"none","text"))
    return(samples)
}

# Function to plot the probability for each cell line of being in a particular state
plot_cell_prob <- function(df, title, HMM_type){
    # i3 or i6 HMM method, need to determine the number of columns on the graph
    df$mag = seq_len(ifelse(HMM_type == "i6", 6, 3))
    long_data <- melt(df, id = "mag")
    long_data$mag <- as.factor(long_data$mag)
    ggplot2::ggplot(long_data, ggplot2::aes_string(x = 'variable', y = 'value', fill = 'mag'))+
        ggplot2::geom_bar(stat="identity", width = 1) +
        ggplot2::coord_flip() +
        ggplot2::theme(
            panel.grid = ggplot2::element_blank(), panel.background = ggplot2::element_blank(),panel.border = ggplot2::element_blank(),
            axis.text=ggplot2::element_text(size=20),
            plot.title = ggplot2::element_text(hjust = 0.5,size = 22),
            #legend.position = "none",
            legend.position="bottom",
            axis.text.x = ggplot2::element_text(size = 16),
            axis.text.y = ggplot2::element_text(size = 16),
            axis.title.x = ggplot2::element_text(size = 18),
            axis.title.y = ggplot2::element_text(size = 18))+
        ggplot2::labs(title = title) +
        #fill = "CNV States") +
        ggplot2::xlab("Cell") +
        ggplot2::ylab("Probability")+
        ggplot2::labs(fill = "States")+
        ggplot2::scale_x_discrete(breaks =seq(1, ncol(df), 9))
}

# Function for total CNV probaility of belonging to each state using THETA prior
cnv_prob <- function(combined_samples) {
    thetas <- combined_samples[,grepl('theta', colnames(combined_samples))]
    #print(paste("Thetas: ", dim(thetas)))
    return(thetas)
}

# Function for each individule cell probabilities, marginalize over the EPSILONS
cell_prob <- function(combined_samples, obj) {
    epsilons <- combined_samples[,grepl('epsilon', colnames(combined_samples))]
    #print(paste("Epsilons: ", dim(epsilons)))
    epsilon_state_frequencies <- apply(as.data.frame(epsilons), 2, function(x) table(factor(x, levels = seq_len(ifelse(getArgs(obj)$HMM_type == "i6", 6, 3)))))
    cell_probs <- epsilon_state_frequencies/colSums(epsilon_state_frequencies)
    return(cell_probs)
}

## Fucntion to Plot the probability of each state for a CNV
plot_cnv_prob <- function(df, title, HMM_type){
    colnames(df) <- seq_len(ifelse(HMM_type == "i6", 6, 3))
    df <- melt(df)
    colnames(df) <- c("row", "State", "Probability")
    states <- as.factor(df$State)
    ggplot2::ggplot(data = df, ggplot2::aes_string(y = 'Probability', x= 'State', fill = 'states')) +
        ggplot2::geom_boxplot()+
        ggplot2::labs(title = title) +
        ggplot2::theme(plot.title = element_text(hjust = 0.5))
}


###############################################
# Main Function to run Bayesion Network Model #
###############################################
#' @title inferCNVBayesNet: Run Bayesian Network Mixture Model To Obtain Posterior Probabilities For HMM Predicted States
#'
#' @description Uses Markov Chain Monte Carlo (MCMC) and Gibbs sampling to estimate the posterior
#' probability of being in one of six Copy Number Variation states (states: 0, 0.5, 1, 1.5, 2, 3) for CNV's identified by
#' inferCNV's HMM. Posterior probabilities are found for the entire CNV cluster and each individual
#' cell line in the CNV.
#'
#' @param file_dir Location of the directory of the inferCNV outputs.
#' @param infercnv_obj InferCNV object.
#' @param HMM_states InferCNV object with HMM states in expression data.
#' @param model_file Path to the BUGS Model file.
#' @param CORES Option to run parallel by specifying the number of cores to be used. (Default: 1)
#' @param out_dir (string) Path to where the output file should be saved to.
#' @param resume_file_token (string) String token that contains some info on settings used to name files.
#' @param postMcmcMethod What actions to take after finishing the MCMC.
#' @param plotingProbs Option for adding plots of Cell and CNV probabilities. (Default: TRUE)
#' @param quietly Option to print descriptions along each step. (Default: TRUE)
#' @param diagnostics Option to plot Diagnostic plots and tables. (Default: FALSE)
#' @param HMM_type The type of HMM that was ran, either 'i3' or 'i6'. Determines how many states were predicted by the HMM.
#' @param k_obs_groups Number of groups in which to break the observations. (default: 1)
#' @param cluster_by_groups If observations are defined according to groups (ie. patients), each group
#'                            of cells will be clustered separately. (default=FALSE, instead will use k_obs_groups setting)
#' @param reassignCNVs (boolean) Given the CNV associated probability of belonging to each possible state, 
#'                            reassign the state assignments made by the HMM to the state that has the highest probability. (default: TRUE)
#' @param no_plot (boolean) Option set by infercnv::run() for producing visualizations.
#' @param useRaster Option to use rasterization when plotting
#'
#' @return Returns a MCMC_inferCNV_obj and posterior probability of being in one of six Copy Number Variation states
#' (states: 0, 0.5, 1, 1.5, 2, 3) for CNV's identified by inferCNV's HMM.
#'
#' @export
#'
#' @examples
#' data(infercnv_data_example)
#' data(infercnv_annots_example)
#' data(infercnv_genes_example)
#' data(infercnv_HMM_states_example)
#'
#' infercnv_object_example <- infercnv::CreateInfercnvObject(raw_counts_matrix=infercnv_data_example, 
#'                                                           gene_order_file=infercnv_genes_example,
#'                                                           annotations_file=infercnv_annots_example,
#'                                                           ref_group_names=c("normal"))
#'           
#' out_dir = tempfile()
#' infercnv_object_example <- infercnv::run(infercnv_object_example,
#'                                          cutoff=1,
#'                                          out_dir=out_dir, 
#'                                          cluster_by_groups=TRUE,
#'                                          analysis_mode="samples",
#'                                          denoise=TRUE,
#'                                          HMM=TRUE,
#'                                          num_threads=2,
#'                                          no_plot=TRUE)
#' mcmc_obj <- infercnv::inferCNVBayesNet(infercnv_obj      = infercnv_object_example,
#'                                        HMM_states        = HMM_states,
#'                                        file_dir          = out_dir,
#'                                        postMcmcMethod    = "removeCNV",
#'                                        out_dir           = out_dir,
#'                                        resume_file_token = "HMMi6.hmm_mode-samples",
#'                                        quietly           = TRUE,
#'                                        CORES             = 2,
#'                                        plotingProbs      = FALSE,
#'                                        diagnostics       = FALSE,
#'                                        HMM_type          = 'i6',
#'                                        k_obs_groups      = 1,
#'                                        cluster_by_groups = FALSE,
#'                                        reassignCNVs      = FALSE,
#'                                        no_plot           = TRUE)
#'                               
inferCNVBayesNet <- function( file_dir,
                              infercnv_obj,
                              HMM_states,
                              out_dir,
                              resume_file_token,
                              model_file        = NULL,
                              CORES             = 1,
                              postMcmcMethod    = NULL,
                              plotingProbs      = TRUE,
                              quietly           = TRUE,
                              diagnostics       = FALSE,
                              HMM_type          = HMM_type,
                              k_obs_groups      = k_obs_groups,
                              cluster_by_groups = cluster_by_groups,
                              reassignCNVs      = TRUE,
                              no_plot           = no_plot,
                              useRaster) {
    
    ################
    # CHECK INPUTS #
    ################
    # if (!file.exists(file_dir)){
    #     error_message <- paste("Cannot find the supplied directory location for the infercnv output.",
    #                            "Please supply the correct path for the output.")
    #     futile.logger::flog.error(error_message)
    #     stop(error_message)
    # }
    if (!is.null(model_file) && !file.exists(model_file)){
        error_message <- paste("Cannot find the model file.",
                               "Please supply the correct path for the model file.")
        futile.logger::flog.error(error_message)
        stop(error_message)
    }
    if (!(CORES == 1)){
        if (is.na(detectCores())){
            futile.logger::flog.warn(paste("Unable to detect number of cores available through parallel:detectCores(), using the provided number", CORES))
        } else if (as.integer(CORES) > detectCores()){
            error_message <- paste("Too many cores previded. The following system has ",detectCores(), " cores.",
                                   "Please select an appropriate amount.")
            futile.logger::flog.error(error_message)
            stop(error_message)
        }
    }
    if(out_dir != "." & !file.exists(out_dir)){
        # create the output directory
        dir.create(file.path(out_dir))
        futile.logger::flog.info(paste("Creating the following Directory: ", out_dir))
    }
    if (is.null(model_file)){
        model_file <- ifelse(HMM_type == "i6", system.file("BUGS_Mixture_Model",package = "infercnv"), system.file("BUGS_Mixture_Model_i3",package = "infercnv"))
    }
    # no_plot will override the plotting options plotingProbs and diagnostics.
    if (no_plot == TRUE) {
        plotingProbs <- FALSE
        diagnostics <- FALSE
    }
        
    args_parsed <- list("file_dir" = file_dir,
                        "model_file" = model_file,
                        "CORES" = CORES,
                        "out_dir"= out_dir,
                        "resume_file_token" = resume_file_token,
                        "plotingProbs" = plotingProbs,
                        "postMcmcMethod"= postMcmcMethod,
                        "quietly" = quietly,
                        "BayesMaxPNormal" = 0,
                        "HMM_type" = HMM_type,
                        "k_obs_groups" = k_obs_groups,
                        "cluster_by_groups" = cluster_by_groups,
                        "reassignCNVs" = reassignCNVs,
                        "diagnostics" = diagnostics)
    #################################
    # LOAD DATA & INITIALIZE OBJECT #
    #################################
    
    ## create the S4 object
    MCMC_inferCNV_obj <- new("MCMC_inferCNV")
    MCMC_inferCNV_obj <- initializeObject(MCMC_inferCNV_obj, args_parsed, infercnv_obj)
    #MCMC_inferCNV_obj <- getStates(MCMC_inferCNV_obj, HMM_states)
    
    #############
    # MEAN & SD #
    #############
    ## Get mean and sd of expression in the predicted cnv areas
    MCMC_inferCNV_obj <- MeanSD(obj = MCMC_inferCNV_obj, 
                                HMM_states = HMM_states, 
                                infercnv_obj = infercnv_obj)
    # check and print the number of genes in each cnv
    if (args_parsed$quietly == FALSE) {
        number_of_genes <- sapply(cellGene(MCMC_inferCNV_obj), function(i) length(i$Genes))
        print(paste("Number of genes for each CNV: ", paste(number_of_genes, sep = " ",collapse = " ")))
        # check the lengths of each cell group
        sapply(cellGene(MCMC_inferCNV_obj), function(x) length(x$Cells))
    }

    ################################
    # Run MCMC Sampling            #
    ################################
    ## Run Gibbs sampling and time the process
    start_time <- Sys.time()
    MCMC_inferCNV_obj <- runMCMC(MCMC_inferCNV_obj, diagnostics)
    end_time <- Sys.time()
    futile.logger::flog.info(paste("Gibbs sampling time: ", difftime(end_time, start_time, units = "min")[[1]], " Minutes"))
    
    ## Save the MCMC.infercnv_object as an RDS
    saveRDS(MCMC_inferCNV_obj, file = file.path(getArgs(MCMC_inferCNV_obj)$out_dir, "MCMC_inferCNV_obj.rds"))
    
    ########
    # Plot #
    ########
    
    # Plot the probability of not being normal state.
    ## Create the title for the plottig the probability of not being normal. 
    if (getArgs(MCMC_inferCNV_obj)$postMcmcMethod == "removeCNV" || getArgs(MCMC_inferCNV_obj)$reassignCNVs == TRUE ){
        title <- sprintf(" (1 - Probabilities of Normal) Before Filtering") 
        output_filename <- "infercnv.NormalProbabilities.PreFiltering"
    } else {
        title <- sprintf(" 1 - Probabilities of Normal ") 
        output_filename <- "infercnv.NormalProbabilities"
    }
    postProbNormal(MCMC_inferCNV_obj,
                   PNormal = NULL,
                   title = title,
                   output_filename = output_filename,
                   useRaster = useRaster)
    
    return(MCMC_inferCNV_obj)
}

#############################################################
# Function to modify CNV's identified base on probabilities #
#############################################################
#' @title filterHighPNormals: Filter the HMM identified CNV's by the CNV's posterior probability
#' of belonging to a normal state.
#'
#' @description The following function will filter the HMM identified CNV's by the CNV's posterior
#' probability of belonging to a normal state identified by the function inferCNVBayesNet(). Will filter
#' CNV's based on a user desired threshold probability. Any CNV with a probability of being normal above
#' the threshold will be removed.
#'
#' @param MCMC_inferCNV_obj MCMC infernCNV object.
#' @param HMM_states InferCNV object with HMM states in expression data.
#' @param BayesMaxPNormal Option to filter CNV or cell lines by some probability threshold.
#' @param useRaster Option to use rasterization when plotting
#'
#' @return Returns a list of (MCMC_inferCNV_obj, HMM_states) With removed CNV's.
#'
#' @export
#' 
#' @examples
#' data(infercnv_mcmc_obj_example)
#' data(infercnv_HMM_states_example)
#' 
#' mcmc_obj_hmm_states_list <- infercnv::filterHighPNormals( MCMC_inferCNV_obj = mcmc_obj, 
#'                                           HMM_states        = HMM_states, 
#'                                           BayesMaxPNormal   = 0.5)
#'

filterHighPNormals <- function( MCMC_inferCNV_obj,
                                HMM_states,
                                BayesMaxPNormal,
                                useRaster) {
    # Add threshold to the MCMC_object
    MCMC_inferCNV_obj <- setBayesMaxPNormal( obj             = MCMC_inferCNV_obj,
                                             BayesMaxPNormal = BayesMaxPNormal )

    ## Either Remove CNV's based on CNV posterier probabilities ("removeCNV")
    ## or remove cell lines based on cell line posterior probabilities ("removeCells")
    if(!(is.null(getArgs(MCMC_inferCNV_obj)$postMcmcMethod))){
        
        if(getArgs(MCMC_inferCNV_obj)$postMcmcMethod == "removeCNV"){
            #MCMC_inferCNV_obj <- removeCNV(MCMC_inferCNV_obj, HMM_states)
            post_removed <- removeCNV(MCMC_inferCNV_obj, HMM_states)
            MCMC_inferCNV_obj <- post_removed[[1]]
            HMM_states <- post_removed[[2]]
        } else {
            MCMC_inferCNV_obj <- removeCells(MCMC_inferCNV_obj, HMM_states)
        }
        
        # Reassign the cnv states based on their probabilities 
        if(getArgs(MCMC_inferCNV_obj)$reassignCNVs == TRUE){
            post_reassign <- reassignCNV(obj        = MCMC_inferCNV_obj, 
                                         HMM_states = HMM_states)
            MCMC_inferCNV_obj <- post_reassign[[1]]
            HMM_states <- post_reassign[[2]]
        }
    }
    
    # Plot the resulting probabilities 
    plotProbabilities(MCMC_inferCNV_obj)
    
    # Plot the Probability of not being normal state
    ## Create the title for the plottig the probability of not being normal 
    if (getArgs(MCMC_inferCNV_obj)$postMcmcMethod == "removeCNV" || getArgs(MCMC_inferCNV_obj)$reassignCNVs == TRUE ){
        title <- sprintf(" (1 - Probabilities of Normal) With Threshold %s", getArgs(MCMC_inferCNV_obj)$BayesMaxPNormal) # MCMC_inferCNV_obj@args$BayesMaxPNormal)
        output_filename <- "infercnv.NormalProbabilities.PostFiltering"
    }
    postProbNormal(MCMC_inferCNV_obj,
                   PNormal = TRUE,
                   title = title,
                   output_filename = output_filename,
                   useRaster = useRaster)
    
    return(list(MCMC_inferCNV_obj, HMM_states))
}

##########################
# Command Line Arguments #
##########################
## Uncomment to use the command line arguments
# if (!is.null(args)){
#     # Set Constants
#     args_parsed <- optparse::parse_args(pargs)
#     file_dir <- args_parsed$file_dir
#     model_file <- args_parsed$model_file
#     CORES <- args_parsed$CORES
#     out_dir <- args_parsed$out_dir
#     inferCNVBayesNet(file_dir,
#                      model_file,
#                      CORES,
#                      out_dir)
# }
