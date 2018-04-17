par(mfrow = c(1,5))
plot_perf(result_BT_10pc, ylim = c(0.1,0.8), main = "Boosting Tree : 10%")
plot_perf(result_RF_10pc, ylim = c(0.1,0.8), main = "Random Forest : 10%")
plot_perf(result_nn_10pc, ylim = c(0.1,0.8), main = "Neural Network : 10%")
plot_perf(result_SVM_10pc, ylim = c(0.1,0.8), main = "SVM : 10%")




metaPred_goodPred = function(rep= 1, outerFold = 1){
    comparePred = result_BT_10pc[[rep]][[outerFold]]$prediction$data[ , c("id", "truth")]
    comparePred = do.call(cbind,
                        args = list(comparePred,
                                    BT_pred = result_BT_10pc[[rep]][[outerFold]]$prediction$data$response,
                                    RF_pred = result_RF_10pc[[rep]][[outerFold]]$prediction$data$response,
                                    NN_pred = result_nn_10pc[[rep]][[outerFold]]$prediction$data$response,
                                    SVM_pred = result_SVM_10pc[[rep]][[outerFold]]$prediction$data$response
                                )
                         )
    
    
    meanPred = c()
    for (i in 1:nrow(comparePred)) {
        meanPred = c(meanPred , names(which.max(table(t(comparePred[i,3:6])))) )
    }
    comparePred = cbind(comparePred, avgPred = meanPred)
    return(sum(as.character(comparePred$truth) == as.character(comparePred$avgPred))/ nrow(comparePred))
}

res = c()
for (i in 1:10) {
    for (j  in 1:3) {
       res = c(res, metaPred_goodPred(rep = i, outerFold = j) )
    }
}
