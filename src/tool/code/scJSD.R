scJSD <- function(data=NULL, classes=NULL){
    print("Inside scJSD func.....")
    if(is.null(data)){ stop("data undefined") }
    if(is.null(classes)){ stop("target class undefined") }   

    col.numbers = which(colnames(data) %in% classes)
    
    ## JSD part
    subpop_names = unique(colnames(data))
    print("printing subpop_name")
    print(subpop_names)
    subpop_names_remaining = setdiff(subpop_names,classes)
    subpop_names = c(classes,subpop_names_remaining) # Placing the target subpopulation as first
    print("printing subpop_names_remaining")
    print(subpop_names_remaining)
    print(subpop_names)
    target_subpopulation_index = 1
    subpopulations_data = lapply(subpop_names, function(subpop_name) {
	data[,grep(paste0(subpop_name,'.*'),colnames(data))]
    })    

    # Remove TFs that are not expressed in any cell in the target
    considered_TFs = rownames(data)        
#    never_expressed_TFs_in_target_subpop = sapply(considered_TFs, function(x) {
#        sum(subpopulations_data[[target_subpopulation_index]][x,]) == 0
#    })
    print("considered_TFs qwerty")
    print(considered_TFs)
    print("subpop_filtered")
    print(subpopulations_data)
    print("subpop_filtered")
    print(subpopulations_data[[target_subpopulation_index]])
    print("subpop_filtered rowMean")
    print(rowMeans(subpopulations_data[[target_subpopulation_index]]))
    print("subpop_filtered rowMean zero check")
    print(rowMeans(subpopulations_data[[target_subpopulation_index]])!=0)
    considered_TFs = considered_TFs[rowMeans(subpopulations_data[[target_subpopulation_index]])!=0]
    print("considered_TFs qwerty3")
    print(!is.na(considered_TFs))
    considered_TFs = considered_TFs[!is.na(considered_TFs)] ##
        
    # Compute background vectors for JSD computation
#    JSD_backround_vectors = as.matrix(data[match(considered_TFs, rownames(data)), -col.numbers, drop=FALSE])
    JSD_backround_vectors = data[match(considered_TFs, rownames(data)), -col.numbers, drop=FALSE]
    # print("Printing bg vectors1")
    # print(col.numbers)
    # print("data hkhlk")
    # print(considered_TFs)
    # print("gfdjdjdjfg")
    # print(rownames(data))
    # print("hkhlk")
    # print(colnames(data))
    # print("ytutyu")
    # print(classes)
    colnames(JSD_backround_vectors) = gsub("^(.+)\\.\\d+$","\\1",colnames(JSD_backround_vectors))
    JSD_backround_vectors = JSD_backround_vectors[,order(colnames(JSD_backround_vectors))]
    # print("Printing bg vectors2")
    # print(JSD_backround_vectors)
    JSD_backround_vectors = t(JSD_backround_vectors)
    # print("Printing bg vectors3")
    # print(JSD_backround_vectors)
    #write.csv(as.matrix(JSD_backround_vectors),"/Users/avani_mahadik/Documents/data_from_code/not_working/JSD_bg_vectors.csv")

    # C++ version
    JSD_values = computeJSDforEachTF(as.matrix(subpopulations_data[[target_subpopulation_index]][considered_TFs,]), as.matrix(JSD_backround_vectors), c(1, rep(0,nrow(JSD_backround_vectors))) )
    # print("printing matrix taken for JSD value calculations")
    # print("matrix 1")
    # print((as.matrix(subpopulations_data[[target_subpopulation_index]][considered_TFs,])))

    # print("matrix 2")
    # print(as.matrix(JSD_backround_vectors))
    
    # print("considered_TFs printing")
    print(col.numbers)
    #write.csv(c(1, rep(0,nrow(JSD_backround_vectors))),"/Users/avani_mahadik/Downloads/scp_pancreatic_islets/ideal.csv")
    #write.csv((subpopulations_data[[target_subpopulation_index]][considered_TFs,]),"/Users/avani_mahadik/Downloads/scp_pancreatic_islets/matrix1.csv")
    #write.csv((subpopulations_data[[target_subpopulation_index]][considered_TFs,]),"/tmp/matrix1.csv")
    #write.csv(considered_TFs,  "/Users/avani_mahadik/Downloads/scp_pancreatic_islets/early_data_transrun/considered_TFs.csv")
    #write.csv(format(data, scientific=FALSE), "/Users/avani_mahadik/Downloads/scp_pancreatic_islets/early_data_transrun/data_mat.csv")
    #write.csv(as.matrix(JSD_backround_vectors),"/Users/avani_mahadik/Downloads/scp_pancreatic_islets/matrix2.csv")
    # write.csv(as.matrix(JSD_backround_vectors),"/Users/avani_mahadik/Documents/data_from_code/not_working/matrix2.csv")



    # R original version
    #D = sapply(1:ncol(b),function(i){i=as.numeric(i); sum(b[,i] - JSD_values[,i])}); sum(abs(D)); which(abs(D) > 1.0e-12)    
#    JSD_values = sapply(considered_TFs, function(tf) {
#        print(tf)
#        #a = computeJSDforEachCell(as.numeric(subpopulations_data[[target_subpopulation_index]][tf,]), c(0,JSD_backround_vectors[,tf]), c(1, rep(0,length(JSD_backround_vectors[,tf]))) )        
#       sapply(subpopulations_data[[target_subpopulation_index]][tf,], function(tf_expr_in_cell) {
#           distr_vector = as.numeric(c(tf_expr_in_cell,JSD_backround_vectors[,tf]))    
#           if(sum(distr_vector) == 0){ return(1) }    
#           distr_vector = distr_vector/sum(distr_vector)    
#           ideal_distr_vector = rep(0,length(distr_vector))
#           ideal_distr_vector[1] = 1    
#           mean_distr_vector = 0.5*(ideal_distr_vector + distr_vector)
#           KL1 = log2(1/mean_distr_vector[1])
#           KL2 = sum( distr_vector * log2(distr_vector / mean_distr_vector),na.rm=TRUE)
#           0.5*(KL1+KL2)    
#       })
#    })

    colnames(JSD_values) = considered_TFs
    rownames(JSD_values) = colnames(subpopulations_data[[target_subpopulation_index]])
    # print("printing JSDvalues")
    # print(JSD_values)
    return(JSD_values)
}

scJSD_score <- function(JSD_values=NULL){
    if(is.null(JSD_values)){ stop("JSD_values undefined") }

    threshold = as.numeric(quantile(as.vector(JSD_values),probs=0.05,type=3))
    # print("**printing threshold**")
    # print(threshold)
    considered_TFs = colnames(JSD_values)
    counts = sapply(considered_TFs, function(tf) {
	sum(JSD_values[,tf] < threshold)#
    })
    print(counts)
    sJSD_expr_score = sort(counts,decreasing=TRUE)
    sJSD_expr_score = sapply(considered_TFs, function(tf) {
        sum(JSD_values[,tf])
    })
    sJSD_expr_score = sort(sJSD_expr_score, decreasing=FALSE) # the smaller the better
    print("Printing sJSD_expr_score")
    print(sJSD_expr_score)
    return(sJSD_expr_score)
}
