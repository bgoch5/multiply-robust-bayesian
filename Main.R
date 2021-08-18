########################
# Main.R #
# Multiply Robust Bayesian Procedures 
# Gochanour et al. (2021+) 


# I. Import Functions, Specify Models & Changes in Functions 
#########################################################
library(dplyr)
source("functions.R")
load("dat.RData")


####################### START MODIFY ################################
# Regression Models
reg1=y~x1+x2+x3+x4+x5
reg2=y~x1*x3+x1*x2+x6+x7+x8

# Propensity Score Models
prop1=r~x1+x2+x3+x4+x5
prop2=r~x1*x3+x1*x2+x6+x7+x8

# Same as above but with (1-r)
prop0_1=(1-r)~x1+x2+x3+x4+x5
prop0_2=(1-r)~x1*x3+x1*x2+x6+x7+x8
####################### END MODIFY ####################################

########## Preparation, Don't Change #############

number=grep("r", colnames(dat))-1

df_reg1 <- subset(dat[,c(1:number,number+2)],dat$r==1)
df_reg1= df_reg1 %>% rename(y=y1)

df_reg0 <- subset(dat[,c(1:number,number+3)],dat$r==0)
df_reg0= df_reg0 %>% rename(y=y0)

df_prop <- dat[,c(1:(number+1))]
#######################################################

# II. Run the Estimators #
###############################################

#######  1. Calibration #######

my_CAL<-fgammaCAL(dat,C=500)
m=my_CAL[1:9]
CrI_lo=my_CAL[10:18]
CrI_hi=my_CAL[19:27]
Res_CAL<-cbind(m,CrI_lo,CrI_hi)
Res_CAL

######## 2. Projection ######

my_PROJ<-fgammaPROJ(dat,C=500)
m=my_PROJ[1:9]
CrI_lo=my_PROJ[10:18]
CrI_hi=my_PROJ[19:27]
Res_PROJ<-cbind(m,CrI_lo,CrI_hi)

# III. Put Together 
###############################################

result=rbind(Res_CAL,Res_PROJ)
rownames(result)= c("CAL_1010",
                    "CAL_0110",
                    "CAL_1001",
                    "CAL_0101",
                    "CAL_1011",
                    "CAL_0111",
                    "CAL_1110",
                    "CAL_1101",
                    "CAL_1111",
                    "PROJ_1010",
                    "PROJ_0110",
                    "PROJ_1001",
                    "PROJ_0101",
                    "PROJ_1011",
                    "PROJ_0111",
                    "PROJ_1110",
                    "PROJ_1101",
                    "PROJ_1111")

finaltable=round(result,2)
finaltable

write.csv(finaltable,"finaltable.csv")
