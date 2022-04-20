
#' @title .i3HMM_get_sd_trend_by_num_cells_fit
#'
#' @description Determines the characteristics for the tumor cell residual intensities, including
#'              fitting the variance in mean intensity as a function of number of cells sampled.
#'
#' @param infercnv_obj infercnv object
#'
#' @param i3_p_val  the p-value to use for defining the position of means for the alternate amp/del distributions.
#'
#' @param plot   boolean, set to TRUE to plot the mean/var fit.
#'
#' @return normal_sd_trend list
#' 
#' @keywords internal
#' @noRd
#'


.i3HMM_get_sd_trend_by_num_cells_fit <- function(infercnv_obj, i3_p_val=0.05, plot=FALSE) {

    
    tumor_samples = infercnv_obj@reference_grouped_cell_indices
    
    tumor_expr_vals <- infercnv_obj@expr.data[,unlist(tumor_samples)]
    
    mu = mean(tumor_expr_vals)
    sigma = sd(tumor_expr_vals)
    nrounds = 100
    sds = c()
    ngenes = nrow(tumor_expr_vals)

    num_tumor_samples = length(tumor_samples)
    flog.info(".i3HMM_get_sd_trend_by_num_cells_fit:: -got ", num_tumor_samples, " samples") 
    
    for (ncells in seq_len(100)) {
        means = c()
        
        for(i in seq_len(nrounds)) {
            ## pick a random gene
            rand.gene = sample(seq_len(ngenes), size=1)
            
            ## pick a random normal cell type
            rand.sample = sample(seq_len(num_tumor_samples), size=1)
            #message("rand.sample: " , rand.sample)
            
            vals = sample(infercnv_obj@expr.data[rand.gene, tumor_samples[[rand.sample]] ], size=ncells, replace=TRUE)
            m_val = mean(vals)
            means = c(means,  m_val)
        }
        sds = c(sds, sd(means))
    }
    
    ## fit linear model
    num_cells = seq_along(sds)
    fit = lm(log(sds) ~ log(num_cells)) #note, hbadger does something similar, but not for the hmm cnv state levels
    
    if (plot) {
        plot(log(num_cells), log(sds), main='log(sd) vs. log(num_cells)')
    }
    
    ## get distribution position according to p_val in a qnorm
    mean_delta = determine_mean_delta_via_Z(sigma, p=i3_p_val)    
    message("mean_delta: ", mean_delta, ", at sigma: ", sigma, ", and pval: ", i3_p_val)
    
    ## do this HBadger style in case that option is to be used.
    KS_delta = get_HoneyBADGER_setGexpDev(gexp.sd=sigma, alpha=i3_p_val)
    message("KS_delta: ", KS_delta, ", at sigma: ", sigma, ", and pval: ", i3_p_val)
    
    
    sd_trend = list(mu=mu,
                    sigma=sigma,
                    fit=fit,
                    mean_delta=mean_delta,
                    KS_delta=KS_delta)
    
    return(sd_trend)
    
}


#' @title .i3HMM_get_HMM
#'
#' @description get the i3 HMM parameterization
#'
#' @param sd_trend  the normal sd trend info list
#'
#' @param num_cells  number of cells in a subcluster
#'
#' @param t  alt state transition probability
#'
#' @param i3_p_val p-value used to define position of mean for amp/del dists (default: 0.05)
#'
#' @param use_KS  boolean : use the KS statistic (HBadger style) for determining the position of the mean for the amp/del dists.
#'
#' @return HMM_info list
#'
#' @noRd

.i3HMM_get_HMM <- function(sd_trend, num_cells, t, i3_p_val=0.05, use_KS) {

    ## Here we do something very similar to HoneyBadger
    ## which is to estimate the mean/var for the CNV states
    ## based on the KS statistic and the expression values
    ## for the normal cells.
    
    
    ##                      states: 0.5,    1,     1 .5
    state_transitions = matrix( c(1-5*t,  t,     t,
                                  t,     1-5*t,  t,
                                  t,      t,     1-5*t),
                               byrow=TRUE,
                               nrow=3)
    
    delta=c(t,     1-5*t,    t) # more likely normal,
    
    mu = sd_trend$mu # normal cell mean
    
    
    #flog.info(sprintf("-.i3HMM_get_HMM, mean_delta=%g, use_KS=%s", mean_delta, use_KS))
    
    if (num_cells == 1) {
        sigma = sd_trend$sigma
        mean_delta = ifelse(use_KS, sd_trend$KS_delta, sd_trend$mean_delta)
    } else {
        ## use the var vs. num cells trend
        sigma <- exp(predict(sd_trend$fit,
                             newdata=data.frame(num_cells=log(num_cells)))[[1]])

        mean_delta = ifelse(use_KS,
                            get_HoneyBADGER_setGexpDev(gexp.sd=sd_trend$sigma, alpha=i3_p_val, k_cells=num_cells), 
                            determine_mean_delta_via_Z(sigma, p=i3_p_val) )
    }
    
    
    state_emission_params = list(mean=c(
                                     mu - mean_delta, # state 0.5 = deletion
                                     mu, # state 1 = neutral
                                     mu + mean_delta), # state 1.5 = amplification
                                 
                                 sd=c(sigma,
                                      sigma,
                                      sigma)  # shared variance as in HB
                                 )
    
    HMM_info = list(state_transitions=state_transitions,
                    delta=delta,
                    state_emission_params=state_emission_params)


    #print(HMM_info)
    
    return(HMM_info)
}


#' @title i3HMM_predict_CNV_via_HMM_on_indiv_cells
#'
#' @description use the i3 HMM for predicting CNV at the level of individual cells
#'
#' @param infercnv_obj infercnv object
#'
#' @param i3_p_val  p-value used to determine mean for amp/del distributions
#'
#' @param sd_trend (optional) by default, computed automatically based on infercnv_obj, i3_p_val
#'
#' @param t alt state transition probability (default: 1e-6)
#'
#' @param use_KS boolean : use the KS test statistic to determine mean for amp/del dist HBadger style (default: TRUE)
#'
#' @return infercnv_obj where infercnv_obj@expr.data contains state assignments.
#'
#' @keywords internal
#' @noRd
#'


i3HMM_predict_CNV_via_HMM_on_indiv_cells  <- function(infercnv_obj,
                                                     i3_p_val=0.05,
                                                     sd_trend=.i3HMM_get_sd_trend_by_num_cells_fit(infercnv_obj, i3_p_val),
                                                     t=1e-6,
                                                     use_KS=TRUE) {
    
    flog.info("predict_CNV_via_HMM_on_indiv_cells()")
    
    chrs = unique(infercnv_obj@gene_order$chr)
    
    expr.data = infercnv_obj@expr.data
    gene_order = infercnv_obj@gene_order
    hmm.data = expr.data
    hmm.data[,] = -1 #init to invalid state

    HMM_info  <- .i3HMM_get_HMM(sd_trend, num_cells = 1, t=t, i3_p_val=i3_p_val, use_KS=use_KS)

    message(HMM_info)
    
    ## run through each chr separately
    lapply(chrs, function(chr) {
        chr_gene_idx = which(gene_order$chr == chr)
        
        ## run through each cell for this chromosome:
        lapply(seq_len(ncol(expr.data)), function(cell_idx) {
            
            gene_expr_vals = as.vector(expr.data[chr_gene_idx,cell_idx])

            hmm <- HiddenMarkov::dthmm(x=gene_expr_vals,
                                       Pi=HMM_info[['state_transitions']],
                                       delta=HMM_info[['delta']],
                                       distn="norm",
                                       pm=HMM_info[['state_emission_params']])
            
            
            hmm_trace <- Viterbi.dthmm.adj(hmm)
            
            hmm.data[chr_gene_idx,cell_idx] <<- hmm_trace
        })
    })
    
    infercnv_obj@expr.data <- hmm.data
    
    return(infercnv_obj)
    
}



#' @title i3HMM_predict_CNV_via_HMM_on_tumor_subclusters
#'
#' @description use the i3 HMM for predicting CNV at the level of tumor subclusters
#'
#' @param infercnv_obj infercnv object
#'
#' @param i3_p_val  p-value used to determine mean for amp/del distributions
#'
#' @param sd_trend (optional) by default, computed automatically based on infercnv_obj, i3_p_val
#'
#' @param t alt state transition probability (default: 1e-6)
#'
#' @param use_KS boolean : use the KS test statistic to determine mean for amp/del dist HBadger style (default: TRUE)
#'
#' @return infercnv_obj where infercnv_obj@expr.data contains state assignments.
#'
#' @keywords internal
#' @noRd
#'

i3HMM_predict_CNV_via_HMM_on_tumor_subclusters  <- function(infercnv_obj,
                                                           i3_p_val=0.05,
                                                           sd_trend=.i3HMM_get_sd_trend_by_num_cells_fit(infercnv_obj, i3_p_val),
                                                           t=1e-6,
                                                           use_KS=TRUE
                                                           ) {
    
    
    flog.info(sprintf("i3HMM_predict_CNV_via_HMM_on_tumor_subclusters(i3_p_val=%g, use_KS=%s)", i3_p_val, use_KS))
    
    if (is.null(infercnv_obj@tumor_subclusters)) {
        flog.warn("No subclusters defined, so instead running on whole samples")
        return(i3HMM_predict_CNV_via_HMM_on_whole_tumor_samples(infercnv_obj, i3_p_val, sd_trend, t, use_KS));
    }
    
    chrs = unique(infercnv_obj@gene_order$chr)
    
    expr.data = infercnv_obj@expr.data
    gene_order = infercnv_obj@gene_order
    hmm.data = expr.data
    hmm.data[,] = -1 #init to invalid state

    tumor_subclusters <- unlist(infercnv_obj@tumor_subclusters[["subclusters"]], recursive=FALSE)
    
    ## add the normals, so they get predictions too:  | already part of tumor_subclusters
    # tumor_subclusters <- c(tumor_subclusters, infercnv_obj@reference_grouped_cell_indices)
    
    ## run through each chr separately
    lapply(chrs, function(chr) {
        chr_gene_idx = which(gene_order$chr == chr)
        
        ## run through each cell for this chromosome:
        lapply(tumor_subclusters, function(tumor_subcluster_cells_idx) {
            
            gene_expr_vals = rowMeans(expr.data[chr_gene_idx,tumor_subcluster_cells_idx,drop=FALSE])
            
            num_cells = length(tumor_subcluster_cells_idx)

            HMM_info  <- .i3HMM_get_HMM(sd_trend, num_cells=num_cells, t=t, i3_p_val=i3_p_val, use_KS=use_KS)
                                                
            hmm <- HiddenMarkov::dthmm(gene_expr_vals,
                                       HMM_info[['state_transitions']],
                                       HMM_info[['delta']],
                                       "norm",
                                       HMM_info[['state_emission_params']])
            
            ## hmm_trace <- HiddenMarkov::Viterbi(hmm)
            hmm_trace <- Viterbi.dthmm.adj(hmm)
            
            hmm.data[chr_gene_idx,tumor_subcluster_cells_idx] <<- hmm_trace
        })
    })
    
    infercnv_obj@expr.data <- hmm.data

    flog.info("-done predicting CNV based on initial tumor subclusters")
        
    return(infercnv_obj)
    
}


#' @title i3HMM_predict_CNV_via_HMM_on_whole_tumor_samples
#'
#' @description use the i3 HMM for predicting CNV at the level of whole tumor samples
#'
#' @param infercnv_obj infercnv object
#'
#' @param i3_p_val  p-value used to determine mean for amp/del distributions
#'
#' @param sd_trend (optional) by default, computed automatically based on infercnv_obj, i3_p_val
#'
#' @param t alt state transition probability (default: 1e-6)
#'
#' @param use_KS boolean : use the KS test statistic to determine mean for amp/del dist HBadger style (default: TRUE)
#'
#' @return infercnv_obj where infercnv_obj@expr.data contains state assignments.
#'
#' @keywords internal
#' @noRd
#'


i3HMM_predict_CNV_via_HMM_on_whole_tumor_samples  <- function(infercnv_obj,
                                                        i3_p_val=0.05,
                                                        sd_trend=.i3HMM_get_sd_trend_by_num_cells_fit(infercnv_obj, i3_p_val),
                                                        t=1e-6,
                                                        use_KS=TRUE
                                                        ) {
    
    
    flog.info("predict_CNV_via_HMM_on_whole_tumor_samples")
    
    
    chrs = unique(infercnv_obj@gene_order$chr)
    
    expr.data = infercnv_obj@expr.data
    gene_order = infercnv_obj@gene_order
    hmm.data = expr.data
    hmm.data[,] = -1 #init to invalid state

    ## add the normals, so they get predictions too:
    tumor_samples <- c(infercnv_obj@observation_grouped_cell_indices, infercnv_obj@reference_grouped_cell_indices)
    
    ## run through each chr separately
    lapply(chrs, function(chr) {
        chr_gene_idx = which(gene_order$chr == chr)
        
        ## run through each cell for this chromosome:
        lapply(tumor_samples, function(tumor_sample_cells_idx) {
            
            gene_expr_vals = rowMeans(expr.data[chr_gene_idx,tumor_sample_cells_idx,drop=FALSE])
            
            num_cells = length(tumor_sample_cells_idx)
            
            HMM_info  <- .i3HMM_get_HMM(sd_trend, num_cells=num_cells, t=t, i3_p_val=i3_p_val, use_KS=use_KS)
            
            hmm <- HiddenMarkov::dthmm(gene_expr_vals,
                                       HMM_info[['state_transitions']],
                                       HMM_info[['delta']],
                                       "norm",
                                       HMM_info[['state_emission_params']])
            
            hmm_trace <- Viterbi.dthmm.adj(hmm)
            
            hmm.data[chr_gene_idx,tumor_sample_cells_idx] <<- hmm_trace
        })
    })
    
    infercnv_obj@expr.data <- hmm.data

    flog.info("-done predicting CNV based on initial tumor subclusters")
        
    return(infercnv_obj)
    
}


#' @title i3HMM_assign_HMM_states_to_proxy_expr_vals
#'
#' @description replace i3 HMM state predictions with their represented CNV levels
#' 
#' @param infercnv_obj infercnv object
#'
#' @return infercnv_obj
#'
#' @keywords internal
#' @noRd
#'


i3HMM_assign_HMM_states_to_proxy_expr_vals <- function(infercnv_obj) {
    
    expr.data = infercnv_obj@expr.data
    
    expr.data[expr.data == 1] <- 0.5
    expr.data[expr.data == 2] <- 1
    expr.data[expr.data == 3] <- 1.5
    
    infercnv_obj@expr.data <- expr.data

    return(infercnv_obj)

}


#' @title determine_mean_delta_via_Z
#'
#' @description determine means for amp/del distributions requiring that they cross the
#'              given distribution based on sigma centered at zero and at the given p value
#'
#' @param sigma standard deviation for a Normal distribution
#'
#' @param p  the p-value at which the distributions should intersect
#'
#' @return delta_for_alt_mean
#'
#' @keywords internal
#' @noRd
#'

determine_mean_delta_via_Z <- function(sigma, p) {
    
    ## want tails of the distribution to minimially overlap at the p-value
    
    delta = abs(qnorm(p=p, mean=0, sd=sigma)) 

    flog.info(sprintf("determine mean delta (sigma: %g, p=%g) -> %g", sigma, p, delta))

    return(delta)
    
}


#' @title get_HoneyBADGER_setGexpDev
#'
#' @description  This method is adapted from HoneyBADGER's setGexpDev method
#'               Essentially, using the KS test to determine where to set the
#'               amp/del means for the distributions.
#'               
#'
#' @param gexp.sd standard deviation for all genes
#'
#' @param alpha the targeted p-value
#'
#' @param k_cells  number of cells to sample from the distribution
#' 
#' @param n_iter number random iterations for sampling from the distribution (default: 100)
#'
#' @param plot boolean, set to True to plot.
#'
#' @return optim.dev
#'
#' @noRd

get_HoneyBADGER_setGexpDev <- function(gexp.sd, alpha, k_cells=2, n_iter=100, plot=FALSE) {

    if (k_cells < 2) {
        flog.warn("get_HoneyBADGER_setGexpDev:: k_cells must be at least 2, setting to 2")
        k_cells = 2
    }
    
    devs <- seq(0, gexp.sd, gexp.sd/10)
    pvs <- unlist(lapply(devs, function(dev) {
        mean(unlist(lapply(seq_len(n_iter), function(i) {
            pv <- ks.test(rnorm(k_cells, 0, gexp.sd), rnorm(k_cells, dev, gexp.sd))
            pv$p.value
        })))
    }))
    if(plot) {
        plot(pvs, devs, xlab="p-value", ylab="deviation", xlim=c(0,1))
    }
    fit <- lm(devs ~ pvs)

    optim.dev <- predict(fit, newdata=data.frame(pvs=alpha))

    #flog.info(sprintf("-get_HoneyBADGER_setGexpDev(sigma=%g, alpha=%g, k_cells=%g) = %g", gexp.sd, alpha, k_cells, optim.dev))
    
    return(optim.dev)
}
