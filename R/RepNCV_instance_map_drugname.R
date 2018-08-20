
# AIM : This function creates some extended version of a Resampling instance by mapping drugnames to a full dataset containing all the dosages
# INPUT : A resampling instance with one dosage for each drugs 
#       + A dataset with several dosages for the same drugs 
# OUTPUT : A new resampling instance, blocked with drug based on the name
RepNCV_instance_map_drugname = function(instance_oneDrug, dataset_oneDosage,  dataset_allDosage, DEBUG = FALSE){

    # THe new resampling ojbject will be like the first one but extended   
    Rep_Nest_CV_instance_allDosage = instance_oneDrug
    
    # For each repetition
    for(i in 1:length(Rep_Nest_CV_instance_allDosage)){
        
        # OUTER FOLD RESAMPLING MAPPING
        # ====== Training sets ======
        # drugNameSet = name of individuals/drugs in the training set based on the ID
        drugNameSet_outer_train = map(get_outer_train_ids(instance_oneDrug, rep_nb = i), function(x){dataset_oneDosage[x, "drugname_typaslab"]})
        # then get all IDs linked to these name in the full dataset
        allDosage_indexes = map(drugNameSet_outer_train, function(x){ which(dataset_allDosage$drugname_typaslab %in% x$drugname_typaslab) } )
        Rep_Nest_CV_instance_allDosage[[i]]$outer$train.inds = allDosage_indexes 
        # ====== Testing sets ======
        # Just use setdiff, what's not in the training set is in the test set
        allDosage_indexes_test = map(allDosage_indexes, function(x){ setdiff(seq(1, nrow(dataset_allDosage)), x ) } )
        Rep_Nest_CV_instance_allDosage[[i]]$outer$test.inds = allDosage_indexes_test
        #Size mapping
        Rep_Nest_CV_instance_allDosage[[i]]$outer$size = nrow(dataset_allDosage)
        
        if(DEBUG){
            print(Rep_Nest_CV_instance_allDosage[[i]]$outer$train.inds)
            a = map(Rep_Nest_CV_instance_allDosage[[i]]$outer$train.inds, function(x){ unique(dataset_allDosage[x , ]$drugname_typaslab) })
            b = map(instance_oneDrug[[i]]$outer$train.inds, function(x){ dataset_oneDosage[x , ]$drugname_typaslab })
            print(map2(.x = a, .y = b, .f = function(x,y) {print(setdiff(y,x))}))
        }
        
        # NEVER MAP INNER IDS DIRECTLY TO THE DATASET, ALWAYS MAP THEM FIRST TO OUTER ID  !!!
        # INNER FOLDS RESAMPLING MAPPING
        for(outer in 1:length(get_outer_train_ids(instance_oneDrug, rep_nb = i))) {
            
            # ====== Training sets ======
            innerID = get_inner_train_ids(NCV_sampling = instance_oneDrug, rep_nb = i, outer_nb = outer)
            outerID = get_outer_train_ids(NCV_sampling = instance_oneDrug, rep_nb = i)[[outer]]
            innerID_to_rawData = map(.x = innerID, .f = function(x){ outerID[x]})
            
            drugNameSet_inner_train = map(innerID_to_rawData, function(x){dataset_oneDosage[x, "drugname_typaslab"]})
            # Here we got names of the drugs that should be in the inner training sets
        
            # Extract a subset of raw dataset based on outer training id and select names of this subset
            dataset_outer_train = dataset_allDosage[get_outer_train_ids(NCV_sampling = Rep_Nest_CV_instance_allDosage, rep_nb = i)[[outer]], ]$drugname_typaslab
            
            allDosage_indexes = map(drugNameSet_inner_train, function(x){ 
                which(dataset_outer_train %in% x$drugname_typaslab) 
            } )
            
            Rep_Nest_CV_instance_allDosage[[i]]$inner[[outer]]$train.inds = allDosage_indexes  
            
            if(DEBUG){
                print(Rep_Nest_CV_instance_allDosage[[i]]$inner[[outer]]$train.inds)
                
                allDosageName_inner = map(.x = get_inner_train_ids(NCV_sampling = Rep_Nest_CV_instance_allDosage, rep_nb = i, outer_nb = outer),
                                          .f = function(x) {
                                              unique(dataset_allDosage[(get_outer_train_ids(Rep_Nest_CV_instance_allDosage, rep_nb = i)[[outer]])[x], ]$drugname_typaslab)
                                          })
                newDrugsName_inner = map(.x = get_inner_train_ids(NCV_sampling = instance_oneDrug, rep_nb = i, outer_nb = outer),
                                         .f = function(x) {
                                             unique(dataset_oneDosage[(get_outer_train_ids(instance_oneDrug, rep_nb = i)[[outer]])[x], ]$drugname_typaslab)
                                         })
                print(map2(.x = allDosageName_inner, .y = newDrugsName_inner, .f = function(x,y) {print(setdiff(y, x))}))
            }
            
            # ====== Testing sets ======
            allDosage_indexes_test = map(allDosage_indexes, function(x){ setdiff(seq(1, length(dataset_outer_train)), x ) } )
            Rep_Nest_CV_instance_allDosage[[i]]$inner[[outer]]$test.inds = allDosage_indexes_test
            Rep_Nest_CV_instance_allDosage[[i]]$inner[[outer]]$size = length(dataset_outer_train)
            
        }
    }
    return(Rep_Nest_CV_instance_allDosage)
}
