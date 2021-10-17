# WtFcpReprocudability.R

# TODO REVISE TO USE NEW DATA SOURCE

library(rtracklayer)
df <- lapply(1:4, function(i)  {
  fn <- sprintf("/n/groups/churchman/rds19/data/S002/features/WT-%d_cliff.gff3", i)
  g <- import(fn)
  g1 <- g[g$type == "cds_segment" & g$seq == 1]
  data.frame(sample = i, 
             gene = as.character(g1$Parent), 
             fcp = width(g1), 
             mean = as.numeric(g1$mean), 
             variance = as.numeric(g1$variance))
})
df <- rbind(df[1][[1]], df[2][[1]], df[3][[1]], df[4][[1]])
df$gene <- factor(df$gene)
df$sample <- factor(df$sample)

x <- xtabs(fcp ~ gene + sample, data = df)
d4 <- apply(x, 1, function(i) ((function(u) u[2] - u[1])(range(i))))
plot.ecdf(d4)
d3 <- lapply(1:4, function(j) apply(x[,-j], 1, function(i) ((function(u) u[2] - u[1])(range(i)))))

par(mfrow = c(2,2), mai = c(0,0,0,0))
for(i in 1:4) {
  plot.ecdf(d3[i][[1]])
  abline(h = c(0.5,0.75, 0.9, 0.95), v = c(10,50, 100), col = "blue")
}

