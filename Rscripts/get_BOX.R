
source('~/Roms_tools/Rscripts/get_roms_data.R')

#Get NO3:
NO3_A1 <- get_box('NO3_roms', avgfile,  c(x1,x2),c(y1,y2))
NO3_A2 <- get_box('NO3_roms', avgfile2, c(x1,x2),c(y1,y2))
NO3_A3 <- get_box('NO3_roms', avgfile3, c(x1,x2),c(y1,y2))
NO3_B1 <- get_box('NO3_roms', avgfile,  c(X1,X2),c(Y1,Y2))
NO3_B2 <- get_box('NO3_roms', avgfile2, c(X1,X2),c(Y1,Y2))
NO3_B3 <- get_box('NO3_roms', avgfile3, c(X1,X2),c(Y1,Y2))

#Get mean size:
L_A1 <- get_box('LNV', avgfile,  c(x1,x2),c(y1,y2))
L_A2 <- get_box('LNV', avgfile2, c(x1,x2),c(y1,y2))
L_A3 <- 3
L_B1 <- get_box('LNV', avgfile,  c(X1,X2),c(Y1,Y2))
L_B2 <- get_box('LNV', avgfile2, c(X1,X2),c(Y1,Y2))
L_B3 <- L_A3

#Get size diversity
Var_A1 <- get_box('VAR', avgfile,  c(x1,x2),c(y1,y2))
Var_A2 <- get_box('VAR', avgfile2, c(x1,x2),c(y1,y2))
Var_B1 <- get_box('VAR', avgfile,  c(X1,X2),c(Y1,Y2))
Var_B2 <- get_box('VAR', avgfile2, c(X1,X2),c(Y1,Y2))

#Get  growth rate at the mean size:
mu_A1 <- get_box('omuNet', bioFile,  c(x1,x2),c(y1,y2))
mu_A2 <- get_box('omuNet', bioFile2, c(x1,x2),c(y1,y2))
mu_B1 <- get_box('omuNet', bioFile,  c(X1,X2),c(Y1,Y2))
mu_B2 <- get_box('omuNet', bioFile2, c(X1,X2),c(Y1,Y2))

#Get d2mudl2
d2mu_A1 <- get_box('od2mudl2', bioFile,  c(x1,x2),c(y1,y2))
d2mu_A2 <- get_box('od2mudl2', bioFile2, c(x1,x2),c(y1,y2))
d2mu_B1 <- get_box('od2mudl2', bioFile,  c(X1,X2),c(Y1,Y2))
d2mu_B2 <- get_box('od2mudl2', bioFile2, c(X1,X2),c(Y1,Y2))

#Calculate muAvg:
muAvg_A1 <- mu_A1$dat + Var_A1$dat*d2mu_A1$dat/2 
muAvg_A2 <- mu_A2$dat + Var_A2$dat*d2mu_A2$dat/2 
muAvg_A3 <- get_box('omuNet', bioFile3, c(x1,x2),c(y1,y2))
muAvg_B1 <- mu_B1$dat + Var_B1$dat*d2mu_B1$dat/2 
muAvg_B2 <- mu_B2$dat + Var_B2$dat*d2mu_B2$dat/2 
muAvg_B3 <- get_box('omuNet', bioFile3, c(X1,X2),c(Y1,Y2))

