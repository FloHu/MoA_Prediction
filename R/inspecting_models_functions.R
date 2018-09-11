plot_reg_path = function(modelInstance){
   lambda_for_pred = modelInstance$learner$par.vals$s
   # Get the index of the lambda used in model building which is the closest to s, the lambda used for testing
   closest_lambda_index = which.min(abs(lambda_for_pred - modelInstance$learner.model$lambda))
   
   coeffs = as.matrix(modelInstance$learner.model$beta[, 1:closest_lambda_index])
   coeffs = coeffs[ rowSums(coeffs)!=0, ] 
   coeffs = as.data.frame(coeffs)
   colnames(coeffs) = seq(1:ncol(coeffs))
   coeffs$feat = rownames(coeffs)
   
   coeffs_plotable = melt(coeffs, id.vars = "feat")
   coeffs_plotable$variable = as.numeric(coeffs_plotable$variable)
   ggplot() + geom_path(data = coeffs_plotable, mapping = aes(x = variable, y = value, group = feat)) +
      geom_hline(yintercept = 0, linetype="dashed") +
      annotate("text", x = ncol(coeffs), y = coeffs[, ncol(coeffs)-1], label = coeffs$feat) +
      theme_bw() + xlim(0,ncol(coeffs)+5)
}
