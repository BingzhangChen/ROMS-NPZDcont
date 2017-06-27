
source('~/Roms_tools/Rscripts/get_roms_data.R')

#Get NO3:
NO3_A1 <- get_boxH('NO3', Hdiv = T, c(x1,x2),c(y1,y2))
NO3_A2 <- get_boxH('NO3', Hdiv = F, c(x1,x2),c(y1,y2))
NO3_B1 <- get_boxH('NO3', Hdiv = T, c(X1,X2),c(Y1,Y2))
NO3_B2 <- get_boxH('NO3', Hdiv = F, c(X1,X2),c(Y1,Y2))

#Get mean size:
L_A1   <- get_boxH('LNV', Hdiv = T, c(x1,x2),c(y1,y2))
L_A2   <- get_boxH('LNV', Hdiv = F, c(x1,x2),c(y1,y2))
L_B1   <- get_boxH('LNV', Hdiv = T, c(X1,X2),c(Y1,Y2))
L_B2   <- get_boxH('LNV', Hdiv = F, c(X1,X2),c(Y1,Y2))

#Get size diversity
Var_A1 <- get_boxH('VAR', Hdiv = T, c(x1,x2),c(y1,y2))
Var_A2 <- get_boxH('VAR', Hdiv = F, c(x1,x2),c(y1,y2))
Var_B1 <- get_boxH('VAR', Hdiv = T, c(X1,X2),c(Y1,Y2))
Var_B2 <- get_boxH('VAR', Hdiv = F, c(X1,X2),c(Y1,Y2))

#Get  growth rate at the mean size:
mu_A1 <- get_boxH('muNet', Hdiv = T, c(x1,x2),c(y1,y2))
mu_A2 <- get_boxH('muNet', Hdiv = F, c(x1,x2),c(y1,y2))
mu_B1 <- get_boxH('muNet', Hdiv = T, c(X1,X2),c(Y1,Y2))
mu_B2 <- get_boxH('muNet', Hdiv = F, c(X1,X2),c(Y1,Y2))

#Get d2mudl2
d2mu_A1 <- get_boxH('d2mudl2', Hdiv = T, c(x1,x2),c(y1,y2))
d2mu_A2 <- get_boxH('d2mudl2', Hdiv = F, c(x1,x2),c(y1,y2))
d2mu_B1 <- get_boxH('d2mudl2', Hdiv = T, c(X1,X2),c(Y1,Y2))
d2mu_B2 <- get_boxH('d2mudl2', Hdiv = F, c(X1,X2),c(Y1,Y2))

#Calculate muAvg:
muAvg_A1 <- mu_A1$dat + Var_A1$dat*d2mu_A1$dat/2 
muAvg_A2 <- mu_A2$dat + Var_A2$dat*d2mu_A2$dat/2 
muAvg_B1 <- mu_B1$dat + Var_B1$dat*d2mu_B1$dat/2 
muAvg_B2 <- mu_B2$dat + Var_B2$dat*d2mu_B2$dat/2 

