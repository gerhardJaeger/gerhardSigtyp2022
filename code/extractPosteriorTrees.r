set.seed(12345)
library(ape)

trees <- c(read.nexus("mrbayes/posterior.run1.t"),
           read.nexus("mrbayes/posterior.run2.t"))


nTrees = min(1000, length(trees))

posteriorSample = sample(trees)[1:nTrees]

write.tree(posteriorSample, "mrbayes/posterior.tree")
