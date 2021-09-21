acc <- function(confu.mat) sum(diag(confu.mat)) / sum(as.vector(confu.mat))
tpr.hit <- function(confu.mat) confu.mat[1,1] / sum(confu.mat[,1])
tpr.miss <- function(confu.mat) confu.mat[2,2] / sum(confu.mat[,2])
fscore <- function(confu.mat) 2 / ((1 / (confu.mat[1,1] / sum(confu.mat[1,]))) + (1 / (confu.mat[1,1] / sum(confu.mat[,1]))))
mcc <- function(confu.mat) (confu.mat[1,1] * confu.mat[2,2] - confu.mat[1,2] * confu.mat[2,1]) / (sqrt(sum(confu.mat[1,])) * sqrt(sum(confu.mat[2,])) * sqrt(sum(confu.mat[,1])) * sqrt(sum(confu.mat[,2])))
auc <- function(roc, h = 0.01) (sum(roc) - (roc[1] + roc[length(roc)]) / 2) * h