instance_creation = function(dataset, printTest = F, nFoldsOuter = 8, nFoldsInner = 8){
    dataset = select(dataset, -drugname_typaslab)
    sampling_method = makeResampleDesc(method = "CV", iters = nFoldsOuter, stratify = TRUE)
    
    # outer instance
    outer_task =  makeClassifTask(data = dataset, target = "process_broad")
    sample_instance_outer =  makeResampleInstance(sampling_method, outer_task)
    
    if(printTest){
        # Testing stratification of outer instance
        for(i in 1:nFoldsOuter){
            cat("Fold", i, "- Train and test sets\n" )
            ind = sample_instance_outer$train.inds[[i]]
            print(table(dataset[ind, "process_broad"]))
            ind = sample_instance_outer$test.inds[[i]]
            print(table(dataset[ind, "process_broad"]))
        }
    }
    
    # defining inner instances
    sampling_method = makeResampleDesc(method = "CV", iters = nFoldsInner, stratify = TRUE)
    sample_instance_inner = list()
    for(i in 1:nFoldsOuter){
        # New task using a subset of the whole dataset defined by outer training set ID
        inner_task = makeClassifTask(data = dataset[ sample_instance_outer$train.inds[[i]] , ], target = "process_broad")
        sample_instance_inner[[i]] = makeResampleInstance(sampling_method, inner_task)
    }
    
    if(printTest){
        # Testing stratification of inner instances
        # BE CAREFUL !!!
        # Indexes of individuals in the inner fold are the indexes of the vector of indexes of the corresponding outer fold 
        for (j in 1:nFoldsOuter ){
            outer_ind = sample_instance_outer$train.inds[[j]]
            cat("Outer fold", j, "\n")
            print(table(dataset[outer_ind, "process_broad"] ))
            cat("Inner folds \n")
            for (i in 1:nFoldsInner){
                ind = c(sample_instance_inner[[j]]$train.inds[[i]], sample_instance_inner[[j]]$test.inds[[i]])
                print(table(dataset[outer_ind[ind], "process_broad"]))
            }
        }
    }
    
    nested_CV_instance = list(outer = sample_instance_outer, inner = sample_instance_inner)
    return(nested_CV_instance)
}
