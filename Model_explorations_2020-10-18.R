library(tidyverse)
library(mlr)

# First need to remember where the trained model is even coming from
# Run the first 140 lines of "Retraining_testset_pred.Rmd"
# Then work with trained_model_loc

trained_model_loc
View(trained_model_loc)

m <- chosen$drug_feature_matrices[[1]]
drugs_v <- m$drugname_typaslab

getLearnerModel(trained_model_loc) # = trained_model_loc$learner.model
getParamSet(getLearnerModel(trained_model_loc))

# this is what we need
getLearnerModel(trained_model_loc$learner.model$next.model)
getOOBPreds(trained_model_loc$learner.model$next.model, task = tsk)

trained_model_loc[["learner.model"]][["next.model"]][["learner.model"]][["oob.times"]]


# get feature importances again:
getLearnerModel(trained_model)
getLearnerModel(trained_model$learner.model$next.model)

# not sure which one of the models I should use???
get_importances <- function (trained_model) {
  mod_imp <- 
    getFeatureImportance(trained_model)$res %>% 
    rename(gene = variable, imp = importance) %>% 
    filter(imp != 0) %>%
    arrange(desc(imp)) %>% 
    mutate(gene = fct_reorder(gene, imp, .desc = TRUE), log2imp = log2(imp)) %>%
    as_tibble()
  
  return(mod_imp)
}  

print(get_importances(trained_model), n = 100)
filter(get_importances(trained_model), gene %in% c("RECG", "DGKA", "MRCB", "NAGA", 
                                                   "YADB", "YIAA", "ASPC"))
get_importances(trained_model$learner.model$next.model)

# apparently no need to worry: the "0" values are simply the variables that were 
# not considered
getFeatureImportance(trained_model)$res

getFeatureImportance(trained_model$learner.model$next.model)$res
getFeatureImportance(trained_model_loc$learner.model$next.model)$res
# one more slot?!
head(trained_model[["learner.model"]][["next.model"]][["learner.model"]][["importance"]])
head(trained_model_loc[["learner.model"]][["next.model"]][["learner.model"]][["importance"]])

m_all <- readRDS("./data/programmatic_output/m_all.rds")
limp <- trained_model_loc[["learner.model"]][["next.model"]][["learner.model"]][["localImportance"]]
colnames(limp) <- paste(m_all$drugname_typaslab, m_all$conc, sep = "_")
limp_bygene <- as_tibble(limp)
limp_bygene <- bind_cols(gene = rownames(limp), limp_bygene)
limp_bygene

limp_bydrug <- as_tibble(t(limp))
limp_bydrug <- bind_cols(m_all[, 1:3], limp_bydrug)
limp_bydrug

# for which drugs is MRCB important?
select(limp_bydrug, drugname_typaslab, conc, process_broad, MRCB) %>% 
  arrange(desc(MRCB)) %>% 
  print(n = 200)

# or vice versa: which genes are driving e.g. MECILLINAM at 0.12 concentration?
# cf. with fingerprint plot: interestingly, mrcB, slt, ddlB, ycfM matter, but 
# NAGA and PAL do not
# at the same time, recG and hupA also contribute - not clear why or how
# hupA: the positive score might help
# recG: perhaps because it's not positive? 
select(limp_bygene, gene, contains("MECILLINAM_0.12")) %>% 
  filter(MECILLINAM_0.12 > 0) %>% 
  arrange(desc(MECILLINAM_0.12)) %>% 
  print(n = 200)

# now what about CEFSULODIN_18?
# again we see all those genes with the same, very low score
# and the only five genes higher than that are:
# mrcB, ycfM, hupA, slt, ybhL
# so in contrast to before, ddlB and recG are now not important - but not clear 
# why 
select(limp_bygene, gene, contains("CEFSULODIN_18")) %>% 
  filter(CEFSULODIN_18 > 0) %>% 
  arrange(desc(CEFSULODIN_18)) %>% 
  print(n = 200)

# and now let's study tunicamcyin, 7.5 - mispredicted
# we see quite different genes now at the top
select(limp_bygene, gene, contains("TUNICAMYCIN_7.5")) %>% 
  filter(TUNICAMYCIN_7.5 > 0) %>% 
  arrange(desc(TUNICAMYCIN_7.5)) %>% 
  print(n = 200)

####################################
######### EXTRACTING TREES #########
####################################

library(randomForest)

forest <- trained_model_loc[["learner.model"]][["next.model"]][["learner.model"]][["forest"]]
View(forest)

my_rf <- trained_model_loc[["learner.model"]][["next.model"]][["learner.model"]]
my_rf
class(my_rf)

### INTERESTING FUNCTIONS in the randomForest package:
# classCenter() to extract prototypes
# importance(my_rf)
# The following two are the same for MeanDecreaseGini but not for the accuracy etc. - ??
head(trained_model_loc[["learner.model"]][["next.model"]][["learner.model"]][["importance"]], 
     n = 2)
head(importance(my_rf), n = 2)
# margin() = ?
# MDSplot() for an MDS plot of the proximity matrix
# outlier() to compute outlying measures based on a proximity matrix (?)
# partialPlot() for a partial dependence plot (?)
# predict() method
# rfcv() - CHECK it!
# study randomForest() documentation to understand what you're getting
# tuneRF() finds the optimal mtry parameter

# Hyperpars were: 200 trees, mtry = 97
getTree(my_rf, k = 1, labelVar = TRUE)

# plotting method - also interesting
plot(my_rf, )
?plot.randomForest

# gives you the size of all of the trees
treesize(my_rf)

# agrees with what we already have in the figure
varImpPlot(my_rf)

# find out which predictor variables are "actually used" (meaning?) in the random forest
used_vars <- varUsed(my_rf)
head(used_vars)

# what it means is how often a variable was used in the forest, i.e. here, for 
# example, 20 times a variable was used once
table(used_vars)

# to get feature indices
my_features <- names(my_rf[["forest"]][["xlevels"]])

# to get the individual variables, we need to be more specific:
# with this call: each column is now one tree, each row a variable:
used_vars_m <- varUsed(my_rf, by.tree = TRUE)
rownames(used_vars_m) <- my_features
dim(used_vars_m)
used_vars_m[1:5, 1:5]
# get the counts over all trees:
var_usg_stats <- 
  apply(used_vars_m, 1, sum) %>% 
  enframe(name = "gene", value = "count") %>% 
  arrange(desc(count))

print(var_usg_stats, n = 50)

# plotting a representative tree:
# see 
# https://stats.stackexchange.com/questions/41443/how-to-actually-plot-a-sample-tree-from-randomforestgettree
# and follow installation instructions

library(randomForest)
library(reprtree)

# so now we get an individual tree:
# compare with console output:
getTree(my_rf, k = 1, labelVar = TRUE)
reprtree:::plot.getTree(my_rf, k = 1, cex = 0.5)

# to get the "representative" tree:
# first need to extract a "representative" tree:
# (not sure why it needs new data)
# again, quite a nuisance to get the actually used data frame? 
my_features
my_dfm <- m_all[, c("drugname_typaslab", "conc", "process_broad", my_features)]

my_reprtree <- 
  reprtree:::ReprTree(rforest = my_rf, newdata = my_dfm[, c("drugname_typaslab", my_features)])

pdf("reprtree.pdf", width = 20, height = 20)
reprtree:::plot.reprtree(my_reprtree)
dev.off()


########################################################
######### IMPLEMENTING MCC (Wang et al 2010) ###########
########################################################

my_tree <- getTree(my_rf, k = 1, labelVar = TRUE)
# to help thinking about it
pdf("tree1.pdf", width = 15, height = 15)
reprtree:::plot.getTree(my_rf, k = 1, cex = 0.5)
dev.off()

my_tree








# 
# 
# ##### **********##### **********##### ********** ##### ********** ##### ********** 
# 
# ###                                       USELESS 
# 
# ##### **********##### **********##### ********** ##### ********** ##### ********** 
# 
# 
# # Hypothesis: our RF model overfits because it does not respect blocking during training
# # IF this is true, then we should see a strong improvement on the test when using 
# # just one concentration as well as a smaller difference between OOB and CV error ...
# # alternatively, we could try to use the caret package
# 
# ### REPEAT FIRST 140 LINES of Retraining_testset_pred with just one concentration per drug ---
# # get the data set
# mc_ext <- readRDS("./data/programmatic_output/mc_ext.rds")
# # CHANGED!
# chosen_one <- filter(mc_ext, fitted_model == "classif.randomForest", 
#                  drug_dosages == "one", !chemical_feats)
# 
# perfs <- 
#   chosen_one$perf_measures[[1]] %>%
#   group_by(cvrep) %>%
#   summarise(mean_mmce = mean(mmce))
# 
# (cv_perf_cplt <- mean(perfs$mean_mmce))
# 
# # our input data set (dfm = drug-feature matrix = same as m_all)
# dfm_one <- chosen_one$drug_feature_matrices[[1]]
# 
# # cvinstance we will use for tuning:
# rin <- chosen_one$resamp_instance[[1]][[1]]
# 
# ## Retraining on the whole training dataset: complete model
# ### Tune parameters
# (lrn <- makeLearner(cl = chosen_one$fitted_model[1], predict.type = "prob"))
# (lrn_wrapped <- makeFilterWrapper(learner = lrn, fw.method = "variance", 
#                                   fw.perc = 0.05))
# 
# (tsk <- make_my_task(dfm = dfm_one, blockvar = "drugname_typaslab"))
# 
# sqrt_p <- floor(sqrt(ncol(dfm_one)))
# rf_grid <- makeParamSet(
#   makeDiscreteParam("ntree", values = c(200, 500, 1000)), 
#   makeDiscreteParam("mtry", values = c(sqrt_p, sqrt_p + 25, sqrt_p + 50, 
#                                        sqrt_p + 75, sqrt_p + 100)), 
#   makeDiscreteParam("fw.perc", values = c(0.25, 0.5, 1))
# )
# 
# (my_file <- "./data/programmatic_output/tuned_params_chosen_one.rds")
# if (!file.exists(my_file)) {
#   tuned_params <- tuneParams(learner = lrn_wrapped, task = tsk, resampling = rin, 
#                              measures = chosen_one$tuning_measure[[1]], par.set = rf_grid, 
#                              control = makeTuneControlGrid())
#   saveRDS(object = tuned_params, file = my_file)
# } else {
#   tuned_params <- readRDS(my_file)
# }
# tuned_params
# 
# # set seed? 
# # original model
# (tuned_lrn <- setHyperPars(lrn_wrapped, par.vals = tuned_params$x))
# 
# # models with proximity matrices, local importances
# (tuned_lrn_prox <- setHyperPars(lrn_wrapped, par.vals = tuned_params$x, 
#                                 proximity = TRUE))
# (tuned_lrn_loc <- setHyperPars(lrn_wrapped, par.vals = tuned_params$x, 
#                                proximity = TRUE, localImp = TRUE))
# 
# train_or_load <- function(tuned_learner, task, modelfile) {
#   if (!file.exists(modelfile)) {
#     message("Modelfile does not exist, training classifier.")
#     trained_model <- mlr::train(learner = tuned_learner, task = task)
#     saveRDS(object = trained_model, file = modelfile)
#     return(trained_model)
#   } else {
#     message("Modelfile already present, loading from file.")
#     return(readRDS(modelfile))
#   }
# }
# 
# trained_model_one <- train_or_load(tuned_learner = tuned_lrn, task = tsk, 
#                                modelfile = "./data/programmatic_output/trained_model_one.rds")
# 
# trained_model_prox <- train_or_load(tuned_learner = tuned_lrn_prox, task = tsk, 
#                                     modelfile = "./data/programmatic_output/trained_model_one_prox.rds")
# 
# trained_model_loc <- train_or_load(tuned_learner = tuned_lrn_loc, task = tsk, 
#                                    modelfile = "./data/programmatic_output/trained_model_one_locimp.rds")
# 
# 
# trained_model_loc
# View(trained_model_loc)
# 
# m <- chosen_one$drug_feature_matrices[[1]]
# drugs_v <- m$drugname_typaslab
# 
# getLearnerModel(trained_model_loc) # = trained_model_loc$learner.model
# getParamSet(getLearnerModel(trained_model_loc))
# 
# # this is what we need
# getLearnerModel(trained_model_loc$learner.model$next.model)
# getOOBPreds(trained_model_loc$learner.model$next.model, task = tsk)
# 
# trained_model_loc[["learner.model"]][["next.model"]][["learner.model"]][["oob.times"]]
# 
# ##### and now it gets interesting, because now we can reevaluate the predictions 
# # from the control set, test set etc. 
# 
# # (starting from line 147 in Retraining_testset_pred.Rmd)
# 
# (split_sets <- readRDS("./data/programmatic_output/split_sets.rds"))
# 
# predict_my_set <- function(trained_model, newdata) {
#   predictions <- predict(trained_model, newdata = as.data.frame(newdata))
#   predictions$data <- as_tibble(
#     cbind(newdata[, c("drugname_typaslab", "conc")], 
#           as.data.frame(predictions))
#   )
#   return(predictions)
# }
# 
# 
# #### Control set
# identical_conds <- split_sets$control_set #%>% 
#   # not sure about the filtering steps
#   #filter(origin == "newsize") # %>% # not sure we should do this filtering step 
#   #filter(paste(drugname_typaslab, conc) %in% 
#    #        paste(split_sets$training_set$drugname_typaslab, split_sets$training_set$conc))
# 
# print(identical_conds, n = 50)
# 
# control_preds <- predict_my_set(trained_model = trained_model_one_loc, 
#                                 newdata = identical_conds)
# performance(control_preds)
# 
# print(control_preds$data, n = 500)
# 
# plot_mcl_probs_heatmap(melt_pred_data_mcl(control_preds$data), mics = mics, 
#                        printplot = FALSE, save = TRUE, file = "./plots/Pred_heatmap_controlset___.pdf")
# 
# 
# #### Test set
# 
# test_preds <- predict_my_set(trained_model = trained_model_one_loc, 
#                              newdata = split_sets$test_set)
# 
# (tst_perf_cplt <- performance(test_preds))
# 
# print(test_preds$data, n = 500)
# 
# plot_mcl_probs_heatmap(melt_pred_data_mcl(test_preds$data), mics = mics, 
#                        printplot = FALSE, save = TRUE, plot_width = 170, plot_height = 125, 
#                        file = "./plots/Pred_heatmap_testset.pdf")
# 
# # usage note:
# # plot_mcl_probs_lines(melt_pred_data_mcl(test_preds$data), save = TRUE, 
# #   printplot = FALSE, file = "./plots/Pred_lines_testset.pdf")
# 
# 
# #### Unknown MoA set
# unknown_preds <- predict_my_set(trained_model = trained_model_one_loc, 
#                                 newdata = split_sets$unknown_moas)
# 
# print(unknown_preds$data, n = 500)
# 
# plot_mcl_probs_heatmap(melt_pred_data_mcl(unknown_preds$data), mics = mics, 
#                        printplot = FALSE, save = TRUE, file = "./plots/Pred_heatmap_unknown.pdf")
# 








