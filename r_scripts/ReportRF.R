library(h2o)
library(tidyverse)
library(gplots)

group = "20190906"
models.directory <- paste0("/n/groups/churchman/rds19/data/S004/h2o/", group, "/")


model.names <- list.files(path = models.directory, pattern = "*.rds")
models <- lapply(model.names, function(infile) {
  readRDS(paste0(models.directory, infile))
})
names(models) <- substr(model.names, start = 1, stop = nchar(model.names) - 4)

results <- as.data.frame(t(sapply(models, function(u) {
  c(fit.gini = h2o.giniCoef(u$fit), fitauc = h2o.auc(u$fit), test.gini = h2o.giniCoef(u$test), testauc = h2o.auc(u$test))
})))

r1 <- results[order(results$test.gini, decreasing = TRUE), ]
print(round(r1, digits = 3))

x <- data.frame(matrix(unlist(strsplit(rownames(r1), "_")), ncol = 3, byrow = TRUE))
names(x) <- c("model", "feature", "window")
r1 <- cbind(r1, x)

x <- reshape2::acast(r1, model ~ window, value.var = "test.gini", fun.aggregate = max)
x <- x[, order(as.integer(colnames(x)))]
x[x <= 0] <- NA
print(x)


heatmap.2(x, 
          Rowv = FALSE,
          Colv = FALSE,
          cellnote = ifelse(is.na(x), "", round(x, 3)),
          margins = c(5, 20),
          col = colorRampPalette(rainbow(4), space = "Lab", bias = 4), 
          dendrogram = "none",
          xlab = "Predictors",
          ylab = "Window size", 
          main = "RF Gini index"
          )

maxdex <- which.max(r1$test.gini)
maxmodel <- rownames(r1)[maxdex]

options(digits = 3)

print(r1)

m_test <- models[[maxmodel]]$test
print(m_test)
m_fit <- models[[maxmodel]]$fit
plot(m_fit, sub = maxmodel)
print(m_fit)
plot(m_test, sub = maxmodel)
h2o.varimp_plot(m_fit, num_of_features = 30)
print(as.data.frame(h2o.varimp(m_fit)))


sessionInfo()
