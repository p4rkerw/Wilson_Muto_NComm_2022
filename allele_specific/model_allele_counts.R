# run interactive
# export LSF_DOCKER_VOLUMES="$HOME:$HOME \
# $STORAGE1/diabneph:$HOME/project \
# $SCRATCH1:$SCRATCH1"
# bsub -Is -G compute-parkerw -n 20 -R 'rusage[mem=128GB]' -q general-interactive -a 'docker(p4rkerw/allele_mod:latest)' /bin/bash

# # run detached
# export LSF_DOCKER_VOLUMES="$HOME:$HOME \
# $STORAGE1/diabneph:$HOME/project \
# $SCRATCH1:$SCRATCH1"
# bsub -G compute-parkerw -J "allele_mods" -n 20 -o "$SCRATCH1/log.allele_mods" -R 'rusage[mem=128GB]' -q general -a 'docker(p4rkerw/allele_mod:latest)' \
# Rscript $SCRATCH1/Wilson_Muto_NComm_2022/allele_specific/model_allele_counts.R

# RUN LOCAL INTERACTIVE R terminal
# SCRATCH1=/mnt/g/scratch
# docker run -it \
# --workdir $HOME \
# -v /mnt/g/diabneph:$HOME/project \
# -v $HOME:$HOME \
# -v $SCRATCH1:$SCRATCH1 \
# -e SCRATCH1="/mnt/g/scratch" \
# p4rkerw/allele_mod:latest

library(data.table)
library(dplyr)
library(tidyr)
library(stringr)
library(Matrix)
library(future.apply)
library(tibble)
library(here)
library(reshape2)
library(broom.mixed)

# read in processed matrices
allele_peak_counts <- fread(here("project","analysis","dkd","allele_specific","allele_counts.csv")) %>%
	as.data.frame()
imprna_allele.mat <- fread(here("project","analysis","dkd","allele_specific","imprna_allele.csv")) %>%
	as.data.frame() %>%
	column_to_rownames(var = "gene")

# enable parallel processing via future package
plan("multicore", workers=20) 
options(future.globals.maxSize = 1280000 * 1024^2) # for 128 Gb RAM
plan()

# create list of unique peaks
peaks <- unique(allele_peak_counts$peak_id)

############
############
# glme model
res.ls <- future_lapply(peaks, function(peak_id){
  peak_mat <- allele_peak_counts[allele_peak_counts$peak_id %in% peak_id,]
  
  # each peak may have multiple sig peak-gene links
  # loop through each combination
  genes <- unique(peak_mat$gene)
  peak_results.ls <- lapply(genes, function(gene) {
	  exp_mat <- imprna_allele.mat[gene,]
	  rownames(exp_mat) <- "imprna"
	  exp_mat <- t(exp_mat) %>% as.data.frame()
	  exp_mat <- rownames_to_column(exp_mat, var="barcode") 
	  mat <- left_join(peak_mat, exp_mat, by="barcode")
	  mat <- mutate(mat, alt_present = ifelse(alt > 0, 1, 0))
	  
	  # normalize the expression interval to 0:100 so a 1 unit increase in gene expression corresponds to 1%
	  mat <- dplyr::mutate(mat, nexp = 100*(imprna - min(imprna))/max(imprna)-min(imprna))
	  
	  # run model
	  fit <- tryCatch({lme4::glmer(y ~ exp + (1|ident), data = data.frame(exp=mat$nexp, y=mat$alt_present, ident=mat$sample), 
	                                      family = binomial(link='logit'))
	  }, error=function(e) NULL)
			  
	  if(is.null(fit)) {
	  	# if model returns ERROR return null
		return(NULL)	
	  } else {
	    res <- tidy(fit, conf.int=TRUE) %>%
		  select(term, estimate, std.error, p.value, conf.low, conf.high) %>%
		  as.data.frame() %>%
		  melt(id.var = "term") %>%
		  dcast(1~variable+term)
	    
            # reorder fixed and random effect columns and format results
	    exp_cols <- colnames(res)[str_ends(colnames(res), "_exp")]
	    ranef_cols <- colnames(res)[str_detect(colnames(res), "sd__")][1] # only take sd of ranef

	    # designate column output order
            output_cols <- unlist(list(exp_cols, ranef_cols))
	    res <- res[,output_cols]
		  
	    res <- cbind(peak=peak_id, gene=gene, res)
	    print(res)
	    return(res)
	  }
   })
  peak_results.df <- bind_rows(peak_results.ls) %>% as.data.frame()
  return(peak_results.df)
})
res <- bind_rows(res.ls) %>% as.data.frame()
write.csv(res, file=here("project","analysis","dkd","allele_specific","glme_results.csv"))

############
############
# glme model with fixed effect for diabetes
res.ls <- future_lapply(peaks, function(peak_id){
  peak_mat <- allele_peak_counts[allele_peak_counts$peak_id %in% peak_id,]
  
  # each peak may have multiple sig peak-gene links
  # loop through each combination
  genes <- unique(peak_mat$gene)
  peak_results.ls <- lapply(genes, function(gene) {
	  exp_mat <- imprna_allele.mat[gene,]
	  rownames(exp_mat) <- "imprna"
	  exp_mat <- t(exp_mat) %>% as.data.frame()
	  exp_mat <- rownames_to_column(exp_mat, var="barcode") 
	  mat <- left_join(peak_mat, exp_mat, by="barcode")
	  mat <- mutate(mat, alt_present = ifelse(alt > 0, 1, 0))
	  
	  # normalize the expression interval to 0:100 so a 1 unit increase in gene expression corresponds to 1%
	  mat <- dplyr::mutate(mat, nexp = 100*(imprna - min(imprna))/max(imprna)-min(imprna))
	  
	  # fit model
	  fit <- tryCatch({lme4::glmer(y ~ exp + diabetes + (1|ident), data = data.frame(exp=mat$nexp, y=mat$alt_present, diabetes=mat$diabetes, ident=mat$sample), 
	                                      family = binomial(link='logit'))
	  }, error=function(e) NULL)
	  
	  if(is.null(fit)) {
	  	# if model returns ERROR return null
		return(NULL)	
	  } else {
	    res <- tidy(fit, conf.int=TRUE) %>%
		  select(term, estimate, std.error, p.value, conf.low, conf.high) %>%
		  as.data.frame() %>%
		  melt(id.var = "term") %>%
		  dcast(1~variable+term)
	    
            # reorder fixed and random effect columns and format results
	    exp_cols <- colnames(res)[str_ends(colnames(res), "_exp")]
	    diabetes_cols <- colnames(res)[str_ends(colnames(res), "_diabetes")]
	    ranef_cols <- colnames(res)[str_detect(colnames(res), "sd__")][1] # only take sd of ranef

	    # designate column output order
            output_cols <- unlist(list(exp_cols, diabetes_cols, ranef_cols))
	    res <- res[,output_cols]
		  
	    res <- cbind(peak=peak_id, gene=gene, res)
	    print(res)
	    return(res)
	  }
   })
  peak_results.df <- bind_rows(peak_results.ls) %>% as.data.frame()
  return(peak_results.df)
})
res <- bind_rows(res.ls) %>% as.data.frame()
write.csv(res, file=here("project","analysis","dkd","allele_specific","glme_dkd_results.csv"))

###############
###############
# glme with diabetes fixed effect and interaction
res.ls <- future_lapply(peaks, function(peak_id){
  peak_mat <- allele_peak_counts[allele_peak_counts$peak_id %in% peak_id,]
  
  # each peak may have multiple sig peak-gene links
  # loop through each combination
  genes <- unique(peak_mat$gene)
  peak_results.ls <- lapply(genes, function(gene) {
	  exp_mat <- imprna_allele.mat[gene,]
	  rownames(exp_mat) <- "imprna"
	  exp_mat <- t(exp_mat) %>% as.data.frame()
	  exp_mat <- rownames_to_column(exp_mat, var="barcode") 
	  mat <- left_join(peak_mat, exp_mat, by="barcode")
	  mat <- mutate(mat, alt_present = ifelse(alt > 0, 1, 0))
	  
	  # normalize the expression interval to 0:100 so a 1 unit increase in gene expression corresponds to 1%
	  mat <- dplyr::mutate(mat, nexp = 100*(imprna - min(imprna))/max(imprna)-min(imprna))
			  
	  # fit model
	  fit <- tryCatch({lme4::glmer(y ~ exp * diabetes + (1|ident), data = data.frame(exp=mat$nexp, y=mat$alt_present, diabetes=mat$diabetes, ident=mat$sample), 
	                                      family = binomial(link='logit'))
	  }, error=function(e) NULL)
	  
	  if(is.null(fit)) {
	  	# if model returns ERROR return null
		return(NULL)	
	  } else {
	    res <- tidy(fit, conf.int=TRUE) %>%
		  select(term, estimate, std.error, p.value, conf.low, conf.high) %>%
		  as.data.frame() %>%
		  melt(id.var = "term") %>%
		  dcast(1~variable+term)
	    
            # reorder fixed and random effect columns and format results
	    exp_cols <- colnames(res)[str_ends(colnames(res), "_exp")]
	    diabetes_cols <- colnames(res)[str_ends(colnames(res), "_diabetes")]
	    interaction_cols <- colnames(res)[str_ends(colnames(res), "_exp:diabetes")]
	    ranef_cols <- colnames(res)[str_detect(colnames(res), "sd__")][1] # only take sd of ranef

	    # designate column output order
            output_cols <- unlist(list(exp_cols, diabetes_cols, interaction_cols, ranef_cols))
	    res <- res[,output_cols]
		  
	    res <- cbind(peak=peak_id, gene=gene, res)
	    print(res)
	    return(res)
	  }
   })
  peak_results.df <- bind_rows(peak_results.ls) %>% as.data.frame()
  return(peak_results.df)
})
res <- bind_rows(res.ls) %>% as.data.frame()
write.csv(res, file=here("project","analysis","dkd","allele_specific","glme_dkdint_results.csv"))
