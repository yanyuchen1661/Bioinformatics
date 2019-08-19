csTest_separate <- function(fitted_model,
                   coef,
                   cell_type = NULL,
                   contrast_matrix = NULL,
                   var_shrinkage = TRUE,
                   sort = T) {
    
    # fitted_model is the output from fitModel().
    # coef is a phenotype name, e.g. "disease",
    #      or a vector of contrast terms, e.g. c("disease", "case", "control").
    # cell_type is a cell type name, e.g. "celltype1", or "neuron".  If cell_type is "JOINT", 
    #      compound effect of coef in all cell types will be tested.
    # contrast_matrix is a matrix (or a vector) to specify contrast, e.g.,
    #      cmat <- matrix(0, 2, 6); cmat[1,3] <- 1: cmat[2,4] <- 1
    #      is to test whether the 3rd parameter and 4th parameter are zero simultaneously
    #      i.e. \beta_{3} = \beta_{4} = 0.
    #      If contrast_matrix is specified, coef and cell_type will be ignored!
    # var_shrinkage is to apply shrinkage on estimated MSE.  Applying shrinkage
    #               helps remove extremely small variance estimation and stablize statistics.
    # sort is a boolean parameter.  The output results will be sorted by p value if sort = TRUE.
    
    N <- fitted_model$N
    K <- ncol(fitted_model$Design_out$Prop)
    W <- fitted_model$Design_out$design_matrix
    beta <- fitted_model$coefs
    beta_var <- fitted_model$coefs_var
    MSE <- fitted_model$MSE
    model_names <- fitted_model$model_names
    G <- length(MSE)
    
    if (is.null(contrast_matrix)) {
        if (tolower(cell_type) == "joint") {
            if (length(coef) == 1) {
                if (length(grep(coef, model_names)) > 0) {
                    cat("Test the joint effect of ", coef, " in all cell types. \n", sep = "")
                    param_vec = grep(coef, model_names)
                    cmatrix = matrix(rep(0, ncol(W) * length(param_vec)), nrow = length(param_vec))
                    for (i in 1:length(param_vec)) {
                        cmatrix[i, param_vec[i]] = 1
                    }
                } else {
                    stop("Coef should be a valid phenotype!")
                }
            } else if (length(coef) == 3) {
                if (length(grep(coef[1], model_names)) >0) {
                    cat("Test the joint effect of ", coef[1], 
                        " level ", coef[2], " vs. level ", coef[3], 
                        " in all cell types. \n", sep = "")
                    param_names = grep(coef[1], model_names, value = T)
                    if (length(grep(paste0(coef[1],coef[3]), param_names)) >0) {
                        control_indx = grep(paste0(coef[1],coef[3]), model_names)
                        if (length(grep(paste0(coef[1],coef[2]), param_names)) >0){
                            case_indx = grep(paste0(coef[1],coef[2]), model_names)
                            cmatrix = matrix(rep(0, ncol(W) * K), nrow = K)
                            for(k in 1:K){
                                cmatrix[k, control_indx[k]] = -1
                                cmatrix[k, case_indx[k]] = 1
                            }
                        } else if (length(grep(paste0(coef[1],coef[2]), param_names)) == 0) {
                            cmatrix = matrix(rep(0, ncol(W) * K), nrow = K)
                            for(k in 1:K) {
                                cmatrix[k, control_indx[k]] = -1
                            }
                        }
                    } else if (length(grep(paste0(coef[1],coef[3]), param_names)) == 0) {
                        if (length(grep(paste0(coef[1],coef[2]), param_names)) == 0) {
                            stop("Contrast levels are not valid!")
                        }
                        case_indx = grep(paste0(coef[1],coef[2]), model_names)
                        cmatrix = matrix(rep(0, ncol(W) * K), nrow = K)
                        for(k in 1:K){
                            cmatrix[k, case_indx[k]] = 1
                        }
                    } else {
                        stop("Coef should be a valid phenotype!")
                    }
                } else {
                    stop("Coef should be a valid phenotype!")
                }
            }
            
        } else if (cell_type %in% colnames(fitted_model$Design_out$Prop)) {
            newK = length(cell_type)
            celltype_names = grep(cell_type, model_names, value = T)
            if (length(coef) == 1) {
                if (length(grep(coef, celltype_names)) > 0){
                    cat("Test the effect of ", coef, " in ", cell_type, ". \n", sep = "")
                    param_vec = match(grep(coef, celltype_names, value = T), model_names)
                    cmatrix = matrix(rep(0, ncol(W) * length(param_vec)), nrow = length(param_vec))
                    for (i in 1:length(param_vec)) {
                        cmatrix[i, param_vec[i]] = 1
                    }
                } else {
                    stop("Coef should be a valid phenotype!")
                }
            } else if (length(coef) == 3) {
                if (length(grep(coef[1], model_names)) >0) {
                    cat("Test the effect of ", coef[1], 
                        " level ", coef[2], " vs. level ", coef[3], 
                        " in ", cell_type,". \n", sep = "")
                    param_names = grep(coef[1], celltype_names, value = T)
                    if (length(grep(paste0(coef[1],coef[3]), param_names)) >0) {
                        control_indx = match(grep(paste0(coef[1],coef[3]), param_names, value = T), 
                                             model_names)
                        if (length(grep(paste0(coef[1],coef[2]), param_names)) >0) {
                            case_indx = match(grep(paste0(coef[1],coef[2]), param_names, value = T), 
                                              model_names)
                            cmatrix = matrix(rep(0, ncol(W) * length(case_indx)), nrow = length(case_indx))
                            for(k in 1:length(case_indx)) {
                                cmatrix[k, control_indx[k]] = -1
                                cmatrix[k, case_indx[k]] = 1
                            }
                        } else if (length(grep(paste0(coef[1],coef[2]), param_names)) == 0) {
                            cmatrix = matrix(rep(0, ncol(W) * length(control_indx)), 
                                             nrow = length(control_indx))
                            for(k in 1:length(control_indx)) {
                                cmatrix[k, control_indx[k]] = -1
                            }
                        }
                    } else if (length(grep(paste0(coef[1],coef[3]), param_names)) == 0) {
                        if (length(grep(paste0(coef[1],coef[2]), param_names)) == 0) {
                            stop("Contrast levels are not valid!")
                        }
                        case_indx = match(grep(paste0(coef[1],coef[2]), param_names, value = T), model_names)
                        cmatrix = matrix(rep(0, ncol(W) * length(case_indx)), nrow = length(case_indx))
                        for (k in 1:length(case_indx)) {
                            cmatrix[k, case_indx] = 1
                        }
                    }
                } else {
                    stop("Coef should be a valid phenotype!")
                }
            } else {
                stop("Coef should be a valid phenotype!")
            }
            
        } else {
            stop("Please specify a valid cell_type value!")
        }
    } else {
        print("contrast_matrix is specified. Coef and cell_type will be ignored.")
        if (is.vector(contrast_matrix)) {
            print("contrast_matrix is a vector.")
            if (length(contrast_matrix) != ncol(W)) {
                stop(paste0("But it should have length ", ncol(W), "!"))
            } else {
                cmatrix <- matrix(contrast_matrix, nrow = 1)
                print(contrast_matrix)
            }
        } else if (is.matrix(contrast_matrix)) {
            print("contrast_matrix is matrix.")
            if (ncol(contrast_matrix) != ncol(W)) {
                stop(paste0("But is should have ", ncol(W), "columns!"))
            } else {
                cmatrix = contrast_matrix
                print(contrast_matrix)
            }
        }
    }
    
    ## apply shrinkage on estimated MSE:
    if (var_shrinkage) {
        MSE_threshold <- quantile(c(MSE), 0.1, na.rm = T)
        MSE[MSE < MSE_threshold] <- MSE_threshold
    }
    
    ## calculate F statistics
    L = cmatrix
    c = rep(0, nrow(L))
    tmp1 <- solve(L %*% solve(t(W) %*% W) %*% t(L))
    inv_sigma <- 1 / MSE
    Lb <-  L %*% beta
    tmp2 <- Lb - matrix(rep(c, ncol(beta)), ncol = ncol(beta))
    a1 <- t(tmp2) %*% tmp1
    
    res_table = data.frame(
        f_statistics = rep(0, G),
        p_value = rep(0, G),
        fdr = rep(0, G)
    )
    res_table$f_statistics <- rowSums(a1 * t(tmp2)) * inv_sigma
    res_table$p_value <- pf(abs(res_table$f_stat), df1 = 1, df2 = N - ncol(W), lower.tail = F)
    res_table$fdr <- p.adjust(res_table$p_value, method = 'fdr')
    
    ## If it is to test one parameter, add beta, beta_var, mu and effect size
    if (nrow(cmatrix) == 1 & sum(cmatrix != 0) == 1 & length(coef) == 1) {
        param_num = which(cmatrix != 0)
        res_table$beta <- beta[param_num, ]
        res_table$beta_var <- beta_var[param_num, ]
        if (param_num > K) {
            res_table$mu <- beta[param_num - K, ]
            res_table$effect_size = res_table$beta / (res_table$mu + res_table$beta / 2)
            res_table <- res_table[, c(4:7, 1:3)]
        }
    }
    
    rownames(res_table) = rownames(fitted_model$Y)
    
    if (sort) {
        res_table = res_table[order(res_table$p_value), ]
    }
    
    return(list(
        res_table = res_table,
        design = fitted_model$Design_out$design
    ))
}
