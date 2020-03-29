
#'---
#'output:
#'  pdf_document:
#'    number_sections: true
#'    toc: true
#'    toc_depth: 2
#'    keep_tex: true
#'title: "nnet models"
#'author: Tim Lucas
#'fontsize: 8pt
#'geometry: margin=0.5in
#'---

#+ defs

name <- 'deep_learning'

#' A fairly standard set of covariates
#' Monthly data.

#+setup, echo = FALSE, cache = FALSE, results = 'hide'

knitr::opts_chunk$set(cache = TRUE, fig.width = 8, fig.height = 5)

set.seed(110120)

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(tidyr)

library(imputeMissings)
library(stringr)
library(parallel)
library(caret)
library(knitr)

source('../../helper_functions/summarise.R')
source('../../helper_functions/extract_year.R')

if(!dir.exists('models')) dir.create('models')


#+ read_pr_data, eval = TRUE

pr <- fread("../../../data/derived/malariaAtlas_pr.csv")



#+ read_back_in

covs_clean <- fread('../../../data/extracted_covs/base.csv')
covs_clean <- covs_clean %>% dplyr::select(-V1)


#+ subset
covs_clean <- covs_clean[pr$year_start >= 2000, ]
pr <- pr %>% filter(year_start >= 2000)


covs_clean <- covs_clean[pr$continent == 'Africa', ]
pr <- pr %>% filter(continent == 'Africa')


#+ trans

covs_clean <-
  covs_clean %>% 
  mutate(accessibility = log1p(accessibility),
         CHIRPS = log1p(CHIRPS),
         VIIRS = log1p(VIIRS))


#+ setup_parallel, cache= FALSE




#'# Base Data
#' ## Random CV
#' 

#+ buikd_caret

keras_custom <- list(label = "Multilayer Perceptron Network with Dropout",
                  library = "keras",
                  loop = NULL,
                  type = c('Regression', "Classification"),
                  parameters = data.frame(
                    parameter = c('epochs', 'size1', 'dropout1', 'size2', 'dropout2', 'size3', 'dropout3', 
                                  "batch_size",
                                  "lr", "rho", "decay", 
                                  "activation1","activation2","activation3",
                                  "activation4"),
                    class = c(rep('numeric', 11), rep("character", 4)),
                    label = c('#Epochs', '#Hidden Units1', 'Dropout Rate1', 
                              '#Hidden Units2', 'Dropout Rate2', 
                              '#Hidden Units3', 'Dropout Rate3', 
                              "Batch Size", "Learning Rate",
                              "Rho", "Learning Rate Decay",
                              "Activation Function1", "Activation Function2", 
                              "Activation Function3", "Output Activation 4")
                  ),
                  grid = function(x, y, len = NULL, search = "grid") {
                    afuncs <- c("sigmoid", "relu", "tanh")
                    if(search == "grid") {
                      out <- expand.grid(
                        size = ((1:len) * 2) - 1, 
                        dropout = seq(0, .7, length = len), 
                        batch_size = floor(nrow(x)/3),
                        lr = 2e-6,
                        rho = .9,
                        decay = 0,
                        activation = "relu"
                      )
                    } else {
                      n <- nrow(x)
                      out <- data.frame(
                        size = sample(2:20, replace = TRUE, size = len),
                        dropout = runif(len, max = .7), 
                        batch_size = floor(n*runif(len, min = .1)),
                        lr = runif(len),
                        rho = runif(len),
                        decay = 10^runif(len, min = -5, 0),
                        activation = sample(
                          afuncs, 
                          size = len, 
                          replace = TRUE
                        )
                      )
                    }
                    out
                  },
                  fit = function(x, y, wts, param, lev, last, classProbs, ...) {
                    require(dplyr)
                    K <- keras::backend()
                    K$clear_session()
                    if(!is.matrix(x)) x <- as.matrix(x)
                    model <- keras::keras_model_sequential()
                    
                    
                    model %>% 
                      keras::layer_dense(
                        units = param$size1, 
                        activation = as.character(param$activation1), 
                        kernel_initializer = keras::initializer_glorot_uniform(),
                        input_shape = ncol(x)
                      ) %>%
                      keras::layer_dropout(rate = param$dropout1,
                                           seed = sample.int(1000, 1)) %>% 
                      keras::layer_dense(
                        units = param$size2, 
                        activation = as.character(param$activation2), 
                        kernel_initializer = keras::initializer_glorot_uniform()
                      ) %>%
                      keras::layer_dropout(rate = param$dropout2,
                                           seed = sample.int(1000, 1))%>% 
                      keras::layer_dense(
                        units = param$size3, 
                        activation = as.character(param$activation3), 
                        kernel_initializer = keras::initializer_glorot_uniform()
                      ) %>%
                      keras::layer_dropout(rate = param$dropout3,
                                           seed = sample.int(1000, 1))
                    
                    
                    
                    if(is.factor(y)) {
                      y <- class2ind(y)
                      model %>% 
                        keras::layer_dense(
                          units = length(lev), 
                          activation = 'softmax'
                        ) %>%
                        keras::compile(
                          loss = "categorical_crossentropy",
                          optimizer = keras::optimizer_rmsprop(
                            lr = param$lr,
                            rho = param$rho,
                            decay = param$decay
                          ),
                          metrics = "accuracy"
                        )
                    } else {
                      model %>% 
                        keras::layer_dense(
                          units = 1, 
                          activation = param$activation4
                        ) %>%
                        keras::compile(
                          loss = "mean_absolute_error",
                          optimizer = keras::optimizer_rmsprop(
                            lr = param$lr,
                            rho = param$rho,
                            decay = param$decay
                          ),
                          metrics = "mean_squared_error"
                        )
                    }
                    model %>% keras::fit(
                      x = x, 
                      y = y,
                      epochs = param$epochs,
                      batch_size = param$batch_size,
                      sample_weight = wts,
                      ...
                    )
                    if(last)
                      model <- keras::serialize_model(model)
                    list(object = model)
                  },
                  predict = function(modelFit, newdata, submodels = NULL) {
                    if(inherits(modelFit$object, "raw"))
                      modelFit$object <- keras::unserialize_model(modelFit$object)
                    if(!is.matrix(newdata)) 
                      newdata <- as.matrix(newdata)
                    out <- predict(modelFit$object, newdata)
                    ## check for model type
                    if(ncol(out) == 1) {
                      out <- out[, 1]
                    } else {
                      out <- modelFit$obsLevels[apply(out, 1, which.max)]
                    }
                    out
                  },
                  prob =  function(modelFit, newdata, submodels = NULL) {
                    if(inherits(modelFit$object, "raw"))
                      modelFit$object <- keras::unserialize_model(modelFit$object)
                    if(!is.matrix(newdata)) 
                      newdata <- as.matrix(newdata)
                    out <- predict(modelFit$object, newdata)
                    colnames(out) <- modelFit$obsLevels
                    as.data.frame(out)
                  },
                  varImp = NULL,
                  tags = c("Neural Network"),
                  sort = function(x) x[order(x$size1, -x$size2),],
                  notes = paste("After `train` completes, the keras model object is serialized",
                                "so that it can be used between R session. When predicting, the", 
                                "code will temporarily unsearalize the object. To make the", 
                                "predictions more efficient, the user might want to use ", 
                                "`keras::unsearlize_model(object$finalModel$object)` in the current", 
                                "R session so that that operation is only done once.",
                                "Also, this model cannot be run in parallel due to",
                                "the nature of how tensorflow does the computations.",
                                
                                "Unlike other packages used by `train`, the `dplyr`",
                                "package is fully loaded when this model is used."),
                  check = function(pkg) {
                    testmod <- try(keras::keras_model_sequential(),
                                   silent = TRUE)
                    if(inherits(testmod, "try-error"))
                      stop("Could not start a sequential model. ",
                           "`tensorflow` might not be installed. ",
                           "See `?install_tensorflow`.", 
                           call. = FALSE)
                    TRUE
                  })


#+ fit_base_random, cache = TRUE, results = 'hide', message = FALSE

n <- length(pr$pf_pr[pr$random_holdout == 0])
len <- 200
afuncs <- c("sigmoid", "relu")

gr <- data.frame(
  batch_size = floor(n*runif(len, min = .05, max = 0.3)),
  lr = runif(len),
  rho = runif(len),
  decay = 10^runif(len, min = -5, 0),
  size1 = sample(10:200, replace = TRUE, size = len),
  dropout1 = runif(len, max = .7), 
  activation1 = sample(
    afuncs, 
    size = len, 
    replace = TRUE
  ),
  size2 = sample(10:200, replace = TRUE, size = len),
  dropout2 = runif(len, max = .7), 
  activation2 = sample(
    afuncs, 
    size = len, 
    replace = TRUE
  ),
  size3 = sample(10:200, replace = TRUE, size = len),
  dropout3 = runif(len, max = .7), 
  activation3 = sample(
    afuncs, 
    size = len, 
    replace = TRUE
  ),
  activation4 = sample(
    c(afuncs, 'linear'), 
    size = len, 
    replace = TRUE
  ),
  epochs = sample(8:50, replace = TRUE, size = len)
  
  
)



m_base_r <- train(y = pr$pf_pr[pr$random_holdout == 0],
                  x = covs_clean[pr$random_holdout == 0, ],
                  method = keras_custom, 
                  weights = pr$examined[pr$random_holdout == 0],
                  tuneGrid = gr,
                  metric = 'MAE',
                  trControl = trainControl(method = 'boot', number = 1, 
                                           search = 'random',
                                           savePredictions = FALSE))

save(m_base_r, file = 'models/base_r.RData')



#+ summary_base_random, cache = FALSE


m_base_r$results %>% 
  dplyr::select(-RMSE, -Rsquared, -RMSESD, -MAESD, -RsquaredSD,
                -activation1, -activation2, -activation3, -activation4) %>% 
  pivot_longer(-MAE) %>% 
  ggplot(aes(x = value, y = MAE)) +
  geom_point() + 
  geom_smooth() + 
  facet_wrap(~ name, scale = 'free') +
  ylim(c(NA, 0.2))

m_base_r$results %>% 
  ggplot(aes(activation1, MAE)) + 
  geom_boxplot() + 
  geom_point() +
  ylim(c(NA, 0.2))
m_base_r$results %>% 
  ggplot(aes(activation2, MAE)) + 
  geom_boxplot() + 
  geom_point() +
  ylim(c(NA, 0.2))
m_base_r$results %>% 
  ggplot(aes(activation3, MAE)) + 
  geom_boxplot() + 
  geom_point() +
  ylim(c(NA, 0.2))

m_base_r$results %>% 
  ggplot(aes(activation4, MAE)) + 
  geom_boxplot() + 
  geom_point() +
  ylim(c(NA, 0.2))


kable(arrange(m_base_r$results[, seq_len(length(m_base_r$bestTune) + 3)], desc(MAE)), digits = 2)

#plot(m_base_r)



#+ predict_base_random, cache = FALSE
pred_base_r <- predict(m_base_r, newdata = covs_clean[pr$random_holdout == 1, ])


summary_base_r <- summarise(pr$pf_pos[pr$random_holdout == 1], 
                            pr$examined[pr$random_holdout == 1],
                            pred_base_r,
                            pr[pr$random_holdout == 1, c('longitude', 'latitude')])

summary <- data.frame(name = paste0('base', name), 
                      covariates = 'base',
                      method = name,
                      cv = 'random',
                      mae = summary_base_r$weighted_mae,
                      correlation = summary_base_r$correlation,
                      time = m_base_r$times$everything[[1]])

write.csv(summary, 'random_base_summary.csv')
summary_base_r$weighted_mae

errors <- data.frame(name = paste0('base', name), 
                     covariates = 'base',
                     method = name,
                     pred = pred_base_r,
                     errors = summary_base_r$errors)

write.csv(errors, 'random_base_errors.csv')




#' next bits

#+ ses

sessionInfo()




