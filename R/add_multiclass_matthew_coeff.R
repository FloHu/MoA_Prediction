
multi_mcc_compute <- function(task, model, pred, feats, extra.args){
    
    conf_mat  = calculateConfusionMatrix(pred)
    conf_mat = conf_mat$result
    classNames =  colnames(conf_mat)
    classNames = classNames[1:length(classNames)-1]
    
    term1 = 0.0
    for(k in classNames){
        for(l in classNames){
            for(m in classNames){
                term1  = term1 + conf_mat[k,k]*conf_mat[l,m] - conf_mat[k,l]*conf_mat[m,k]
            }
        }
    }
    
    
    term2 = 0.0
    for(k in classNames){
        term21 = 0.0
        term22 = 0.0
        for(l in classNames){
            term21 = term21 + conf_mat[k,l]
        }
        for (kp in setdiff(classNames, k)){
            for(lp in classNames){
                term22  = term22 + conf_mat[kp,lp]
            }
        }
        term2 = term2 + term21*term22
    }
    
    term3 = 0.0
    for(k in classNames){
        term31 = 0.0
        term32 = 0.0
        for(l in classNames){
            term31 = term31 + conf_mat[l,k]
        }
        for (kp in setdiff(classNames, k)){
            for(lp in classNames){
                term32  = term32 + conf_mat[lp,kp]
            }
        }
        term3 = term3 + term31*term32
    }
    
    
    mcc = term1 / (sqrt(term2)*sqrt(term3))
    return(mcc)
}



multiclass_mcc <- makeMeasure(
    id = "multiclass_mcc", name = "Multiclass Matthew Correlation Coefficient",
    properties = c("classif", "classif.multi"),
    minimize = T,
    best = 1, worst = -1,
    fun = multi_mcc_compute
)


performance(res$pred, multiclass_mcc)


