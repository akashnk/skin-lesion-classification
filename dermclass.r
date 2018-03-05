derm_data = read.table("/media/Omega/Development/R/New/Dermatology/dermatology.data.txt" , header = FALSE , sep = ",")
derm_data <- derm_data[,c(1:33,35)]

colnames(derm_data) <- c("erythema",
                         "scaling",
                         "definite borders",
                         "itching",
                         "koebner phenomenon",
                         "polygonal papules",
                         "follicular papules",
                         "oral mucosal involvement",
                         "knee and elbow involvement",
                         "scalp involvement",
                         "family history, (0 or 1)",
                         "melanin incontinence",
                         "eosinophils in the infiltrate",
                         "PNL infiltrate",
                         "fibrosis of the papillary dermis",
                         "exocytosis",
                         "acanthosis",
                         "hyperkeratosis",
                         "parakeratosis",
                         "clubbing of the rete ridges",
                         "elongation of the rete ridges",
                         "thinning of the suprapapillary epidermis",
                         "spongiform pustule",
                         "munro microabcess",
                         "focal hypergranulosis",
                         "disappearance of the granular layer",
                         "vacuolisation and damage of basal layer",
                         "spongiosis",
                         "saw-tooth appearance of retes",
                         "follicular horn plug",
                         "perifollicular parakeratosis",
                         "inflammatory monoluclear inflitrate",
                         "band-like infiltrate",
                         "classes"
                         )

derm_data$classes <- ifelse(derm_data$classes == "1", "psoriasis",
                          ifelse(derm_data$classes == "2", "seboreic_dermatitis",
                                 ifelse(derm_data$classes == "3", "lichen_planus",
                                        ifelse(derm_data$classes == "4", "pityriasis_rosea",
                                               ifelse(derm_data$classes == "5", "chronic_dermatitis", "pityriasis_rubra_pilaris")))))
                                                      

#derm_data[derm_data == "?"] <- NA

library(mice)

derm_data[,1:33] <- apply(derm_data[, 1:33], 6, function(x) as.numeric(as.character(x)))
dataset_impute <- mice(derm_data[, 1:33],  print = FALSE)
derm_data <- cbind(derm_data[, 34, drop = FALSE], mice::complete(dataset_impute, 1))

derm_data$classes <- as.factor(derm_data$classes)

library(ggplot2)
p <- ggplot(derm_data, aes(x = classes, fill = classes))
p + geom_bar()

ggplot(derm_data, aes(x = erythema)) +
  geom_histogram(bins = 3)

library(pcaGoPromoter)
library(ellipse)
library(lattice)

# perform pca and extract scores
#suppressWarnings(as.numeric(c("1", "2", "X")))

pcaOutput <- pca(t(derm_data[, -1]), printDropped = FALSE, scale = TRUE, center = TRUE)
pcaOutput2 <- as.data.frame(pcaOutput$scores)

# define groups for plotting
pcaOutput2$groups <- derm_data$classes

centroids <- aggregate(cbind(PC1, PC2) ~ groups, pcaOutput2, mean)

conf.rgn  <- do.call(rbind, lapply(unique(pcaOutput2$groups), function(t)
  data.frame(groups = as.character(t),
             ellipse(cov(pcaOutput2[pcaOutput2$groups == t, 1:2]),
                     centre = as.matrix(centroids[centroids$groups == t, 2:3]),
                     level = 0.95),
             stringsAsFactors = FALSE)))

ggplot(data = pcaOutput2, aes(x = PC1, y = PC2, group = groups, color = groups)) +
  geom_polygon(data = conf.rgn, aes(fill = groups), alpha = 0.2) +
  geom_point(size = 2, alpha = 0.6) +
  scale_color_brewer(palette = "Set1") +
  labs(color = "",
       fill = "",
       x = paste0("PC1: ", round(pcaOutput$pov[1], digits = 2) * 100, "% variance"),
       y = paste0("PC2: ", round(pcaOutput$pov[2], digits = 2) * 100, "% variance"))


library(tidyr)

gather(derm_data, x, y, erythema:'knee and elbow involvement') %>%
  ggplot(aes(x = y, color = classes, fill = classes)) +
  geom_density(alpha = 0.3) +
  facet_wrap( ~ x, scales = "free", ncol = 3)

gather(derm_data, x, y, 'scalp involvement':'hyperkeratosis') %>%
  ggplot(aes(x = y, color = classes, fill = classes)) +
  geom_density(alpha = 0.3) +
  facet_wrap( ~ x, scales = "free", ncol = 3)

gather(derm_data, x, y, 'parakeratosis':'band-like infiltrate') %>%
  ggplot(aes(x = y, color = classes, fill = classes)) +
  geom_density(alpha = 0.3) +
  facet_wrap( ~ x, scales = "free", ncol = 4)

# configure multicore
library(doParallel)
cl <- makeCluster(detectCores())
registerDoParallel(cl)

library(caret)

set.seed(42)
index <- createDataPartition(derm_data$classes, p = 0.7, list = FALSE)
train_data <- derm_data[index, ]
test_data  <- derm_data[-index, ]
set.seed(42)

library(dplyr)

rbind(data.frame(group = "train", train_data),
      data.frame(group = "test", test_data)) %>%
  gather(x, y, erythema:'parakeratosis') %>%
  ggplot(aes(x = y, color = group, fill = group)) +
  geom_density(alpha = 0.3) +
  facet_wrap( ~ x, scales = "free", ncol = 3)

model_glm <- caret::train(erythema ~ .,
                          data = train_data,
                          method = "glm",
                          preProcess = c("scale", "center"),
                          na.action=na.exclude,
                          trControl = trainControl(method = "repeatedcv", 
                                                   number = 10, 
                                                   repeats = 10, 
                                                   savePredictions = TRUE, 
                                                   verboseIter = FALSE))

model_glm



set.seed(42)
model_derm <- caret::train(classes ~ .,
                         data = train_data,
                         method = "derm",
                         preProcess = c("scale", "center"),
                         na.action=na.exclude,
                         trControl = trainControl(method = "repeatedcv", 
                                                  number = 10, 
                                                  repeats = 10, 
                                                  savePredictions = TRUE, 
                                                  verboseIter = FALSE))

model_derm$finalModel$confusion
imp <- model_derm$finalModel$importance
imp[order(imp, decreasing = TRUE), ]

# estimate variable importance
importance <- varImp(model_derm, scale = TRUE)
plot(importance)

confusionMatrix(predict(model_derm, test_data), test_data$classes)

results <- data.frame(actual = test_data$classes,
                      predict(model_derm, test_data, type = "prob"))

results$prediction <- ifelse(results$pityriasis_rosea > 1/6, "pityriasis_rosea",
                             ifelse(results$chronic_dermatitis > 1/6, "chronic_dermatitis",
                                    ifelse(results$lichen_planus > 1/6, "lichen_planus",
                                           ifelse(results$pityriasis_rubra_pilaris > 1/6, "pityriasis_rubra_pilaris",
                                                  ifelse(results$psoriasis > 1/6, "psoriasis","seboreic_dermatitis")))))
                                                         

results$correct <- ifelse(results$actual == results$prediction, TRUE, FALSE)

ggplot(results, aes(x = prediction, fill = correct)) +
  geom_bar(position = "dodge")

ggplot(results, aes(x = prediction, y = psoriasis, color = correct, shape = correct)) +
  geom_jitter(size = 3, alpha = 0.6)

ggplot(results, aes(x = prediction, y = chronic_dermatitis, color = correct, shape = correct)) +
  geom_jitter(size = 3, alpha = 0.6)

ggplot(results, aes(x = prediction, y = lichen_planus, color = correct, shape = correct)) +
  geom_jitter(size = 3, alpha = 0.6)

ggplot(results, aes(x = prediction, y = pityriasis_rubra_pilaris , color = correct, shape = correct)) +
  geom_jitter(size = 3, alpha = 0.6)

ggplot(results, aes(x = prediction, y = pityriasis_rosea, color = correct, shape = correct)) +
  geom_jitter(size = 3, alpha = 0.6)

ggplot(results, aes(x = prediction, y = seboreic_dermatitis, color = correct, shape = correct)) +
  geom_jitter(size = 3, alpha = 0.6)


