
get_innerCV = function(repeat_ind, data = Rep_Nest_CV_instance){
    return(data[[repeat_ind]]$inner)
}

get_outerCV = function(repeat_ind, data = Rep_Nest_CV_instance){
    return(data[[repeat_ind]]$outer)
}
