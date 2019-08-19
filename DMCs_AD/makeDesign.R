makeDesign <- function(design, Prop) {

    # design is a N by P matrix, with rows as samples and columns as phenotypes (e.g. age, gender, disease, etc.)
    # Prop is a N by K matrix, with rows as samples and columns as cell types

    K = ncol(Prop)

    if (is.null(colnames(Prop))) {
        colnames(Prop) <- paste0("celltype", 1:K)
    }
    dd <- cbind(Prop, design)

    ## make a formula
    formul <- paste("~", paste(colnames(dd)[1:K], collapse="+"))
    for(i in (K+1):ncol(dd)) {
        ## interaction terms for this factor
        intterms <- paste( paste(colnames(dd)[1:K], colnames(dd)[i], sep=":"), collapse="+")
        formul <- paste(formul, intterms, sep="+")
    }
    design_matrix <- model.matrix(as.formula(formul), dd)[,-1]

    ## remove columns that are all zero
    design_matrix <- design_matrix[,!colSums(design_matrix==0) == nrow(Prop)]

    return(list(design_matrix = design_matrix,
                Prop = Prop,
                design = design,
                all_coefs = colnames(design),
                all_cell_types = colnames(Prop),
                formula = formul))
}
