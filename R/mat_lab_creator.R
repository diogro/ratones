source('./R/read_ratones.R')
library(R.matlab)
library(MasterBayes)
 
ped <- unique(read.csv("./data/pedigree.csv"))
#fixing pedigree so ID are not repeated
count_ID <- table(ped[,1])
ped <- ped[!(ped[,1] %in% names(which(count_ID > 1)) & is.na(ped[,2])),]
ped = orderPed(ped)

if(!require(pedRSP)) install_github("jraffa/pedRSP");
library(pedRSP)
A = solve(Ainv[pos, pos])
rsp = computeRSPs(kmat = A, thre = 0.5)

Ainv = as.matrix(inverseA(ped)$Ainv)
pos = aaply(full_data$ID, 1, function(x) which(x == ped$ID))

Z = matrix(0, nrow(ped), nrow(Y))
for(col in 1:nrow(Y))
  Z[pos[col], col] = 1

Y = as.matrix(dplyr::select(full_data, IS_PM:BA_OPI))
X = model.matrix(lm(Y ~ full_data$SEX + full_data$LIN + full_data$AGE))

writeMat("./Gmatlab/ratones/setup.mat", A = solve(Ainv), X = t(X), Y = Y, Z_1 = Z)

######
# Control
#####

strain = "h'"
for(strain in unique(full_data$line)){
    current_data = full_data %>% filter(line == strain)
    Y = as.matrix(current_data %>% dplyr::select(IS_PM:BA_OPI))
    pos = aaply(current_data$ID, 1, function(x) which(x == ped$ID))
    Z = matrix(0, nrow(ped), nrow(Y))
    for(col in 1:nrow(Y))
        Z[pos[col], col] = 1
    X = model.matrix(lm(Y ~ current_data$SEX + current_data$AGE))
    writeMat(paste0("./Gmatlab/", strain, "/setup.mat"), A = solve(Ainv), X = t(X), Y = Y, Z_1 = Z)
}



library(ggplot2)
library(reshape2)
if(!require(viridis)) install.packages("viridis")
library(viridis)
plotMatrix <- function (corMat, file = NULL) {
  diag(corMat) <- NA
  n_traits = nrow(corMat) 
  myPalette <- viridis(50)
  ## Se quiser uma paleta All American, use essa linha em vez da anterior
  #myPalette <- colorRampPalette(c("blue", "white", "red"))(n = 50)
  m.rs = melt(corMat) 
  m.rs$Var1 <- factor(m.rs$Var1, levels = m.rs$Var1[n_traits:1])
  m.rs.position = m.rs
  m.rs.position$Var1 <- as.numeric(m.rs.position$Var1)
  m.rs.position$Var2 <- as.numeric(m.rs.position$Var2)
  m.rs.position$value= round(m.rs.position$value, 2)
  m.rs.position$value[is.na(m.rs.position$value)] <- levels(m.rs$Var1)[n_traits:1]
  p <- 
    ggplot (m.rs) +
    geom_tile(aes(x = Var2, y = Var1, fill = value)) +
    scale_fill_gradientn(name = '', colours = myPalette) +
    labs(x = NULL, y = NULL) + 
    geom_text(data = m.rs.position, aes(x = Var2, y = Var1, label = value)) + 
    theme_bw()
  if(!is.null(file)) cowplot::save_plot(plot = p, file)
  return(p)
}
plotMatrix(cov2cor(P))

folder = "t"
readMatLab <- function(folder){
    x = readMat(paste0("./Gmatlab/", folder, "/Posterior_mean.mat"))
    G = x$posterior.mean[,,1]$G
    P = x$posterior.mean[,,1]$P
    Gs = x$posterior.mean[,,1]$G
    Ps = x$posterior.mean[,,1]$P
    list(P = P, Ps = Ps, G = G, Gs = Gs)
    }
readMatLab("t")

x = readMat("./Gmatlab/ratones/Posterior_mean.mat")
G = x$posterior.mean[,,1]$G
P = x$posterior.mean[,,1]$P
MatrixCompare(G, Wmat)
MatrixCompare(P, Wmat)
MatrixCompare(P, G)
MatrixCompare(Calomys.Pacote$G, Calomys.Pacote$P)

RandomSkewers(llply(r_models, function(x) x$MAP), G)
KrzCor(llply(r_models, function(x) x$MAP), G)
P
G
Wmat
diag(G)/diag(P)
