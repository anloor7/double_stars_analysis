

########## DATA ANALYSIS ##########


# Filtering DMSA data. First we retain only double systems. Then we remove the systems with 'A' or 'B' in column 
# Qual and the systems whose angular separation is below 3 arcsec. In addition, we add deltaHP column


dmsa_data_previous <- read.csv("dmsa_data.csv", 
                               header=TRUE)
head(dmsa_data_previous)
nrow(dmsa_data_previous)
ncol(dmsa_data_previous)
dmsa_data <- dmsa_data_previous[,-9] # We remove the last column, which has not meaning 
head(dmsa_data)


# Filtering 

dmsa_data$Ncomp <= 2
dmsa_doubles_data <- subset(dmsa_data, Ncomp <= 2) # Retaining only the double systems 
head(dmsa_doubles_data, 4)
nrow(dmsa_doubles_data)
reliable_dmsa_doubles_data <- subset(dmsa_doubles_data, Qual=='A' | Qual=='B') # Filtering by Qual
auxiliary_rho <- seq(1, nrow(reliable_dmsa_doubles_data), by = 2)
reliable_dmsa_doubles_data$rho[auxiliary_rho] <- reliable_dmsa_doubles_data$rho[auxiliary_rho + 1] # Completing NAs in rho 
filtered_dmsa_doubles_data <- subset(reliable_dmsa_doubles_data, rho < 3) # Filtering by rho 
nrow(filtered_dmsa_doubles_data)

Hpdoubles <- filtered_dmsa_doubles_data$Hp
length(Hpdoubles)
head(filtered_dmsa_doubles_data, 2)

auxiliary_Hp <- numeric(length(Hpdoubles))
for (i in seq(1, length(auxiliary_Hp), by = 2)) {
  auxiliary_Hp[i] <- Hpdoubles[i+1]-Hpdoubles[i]
  auxiliary_Hp[i+1] <- auxiliary_Hp[i]
}

filtered_dmsa_doubles_data$deltaHp <- auxiliary_Hp # Creating the column deltaHp 
reduced_filtered_dmsa_doubles_data <- filtered_dmsa_doubles_data[seq(1, 
                                                                     nrow(filtered_dmsa_doubles_data), by = 2),]# Dataframe with only one component per system (the A component)


nrow(filtered_dmsa_doubles_data)
nrow(reduced_filtered_dmsa_doubles_data) # We have 7316 doubles in the sample 

# Obtaining the values of astrometric paremeters in Hipparcos catalogue 

hip <- reduced_filtered_dmsa_doubles_data$HIP
auxiliary_hip <- numeric(length(hip))
for (i in 1:length(auxiliary_hip)) {
  auxiliary_hip[i] <- (paste0(as.numeric(hip[i]), ','))
}


writeLines(auxiliary_hip, sep = '') # We copy these hip numbers from the console and make the query for astrometric quantities in 
# the Hipparcos catalogue 


# Loading Hipparcos data and merging with dmsa_data 

hipparcos_data <- read.csv('hipparcos_data.csv', header = T)
head(hipparcos_data)
merge_hipparcos_dmsa <- merge(reduced_filtered_dmsa_doubles_data, hipparcos_data)
head(merge_hipparcos_dmsa)

library(FITSio)
xm <- readFITS(file = "HGCA_vDR2_corrected.fits", hdu = 1, maxLines = 500000)
head(xm$col)
xm$colNames # We want columns 1, 2, 8, 9, 10, 11, 12, 13, 14

auxiliary_hip <- which(xm$col[[1]] %in% hip) # 5701 stars out of 7316 in the XM 
HIP <- xm$col[[1]][auxiliary_hip]
gaia_id <- xm$col[[2]][which(xm$col[[1]] %in% hip)]
gaia_parallax <-xm$col[[8]][which(xm$col[[1]] %in% hip)]
gaia_e_parallax <- xm$col[[9]][which(xm$col[[1]] %in% hip)]
gaia_pmra <- xm$col[[10]][which(xm$col[[1]] %in% hip)]
gaia_pmde <- xm$col[[11]][which(xm$col[[1]] %in% hip)]
gaia_pmra_e <- xm$col[[12]][which(xm$col[[1]] %in% hip)]
gaia_pmde_e <- xm$col[[13]][which(xm$col[[1]] %in% hip)]
gaia_pmra_pmde <- xm$col[[14]][which(xm$col[[1]] %in% hip)]

gdr2_df <- data.frame(HIP, gaia_id, gaia_parallax, gaia_e_parallax, gaia_pmra,  # Data about GDR2 
                      gaia_pmde, gaia_pmra_e, gaia_pmde_e, gaia_pmra_pmde)
head(gdr2_df)
nrow(gdr2_df)

hipparcos_gdr2 <- merge(merge_hipparcos_dmsa, gdr2_df, by = 'HIP')
head(hipparcos_gdr2)
nrow(hipparcos_gdr2)

# Parallax analysis 

omega_h <- hipparcos_gdr2$Plx
omega_g <- hipparcos_gdr2$gaia_parallax
sigma_h <- hipparcos_gdr2$e_Plx
sigma_g <- hipparcos_gdr2$gaia_e_parallax
var_h <- sigma_h^2
var_g <- sigma_g^2
vector_muestral <- (omega_g-omega_h)^2/(var_g+var_h)
ks.test(vector_muestral, "pchisq", 1)

library(ggplot2)
library(car)
library(tidyverse)

qqPlot(vector_muestral, dist='chisq', df=1, ylim=c(0, 20), col='orange', xlab='Theoretical Quantiles', ylab='Sample Quantiles',
       main='Parallaxes', grid=F, envelope = 0.95)
df_vector_muestral <- as.data.frame(vector_muestral)
ggplot(data.frame(df_vector_muestral), aes(sample = df_vector_muestral$vector_muestral)) + 
  stat_qq(distribution = stats::qchisq, dparams = list(df = 1), col = 'orange') + 
  stat_qq_line(distribution = stats::qchisq, dparams = list(df = 1), col = 'blue')

ranks <- rank(vector_muestral, ties.method="max")
F_plot <- (ranks-0.5)/length(vector_muestral) 
F0inv <- qchisq(F_plot, df = 2)
df_plot <- data_frame(vector_muestral, F0inv)
p <- ggplot(df_plot, aes(F0inv, vector_muestral))
p + geom_point(col = 'blue') + geom_abline(intercept = 0, slope = 1, col = 'orange') +
  xlab('Theoretical quantiles') + ylab('Sample quantiles') + labs(title = 'Parallaxes')

# Proper motions analysis

pmra_g_vector <- hipparcos_gdr2$gaia_pmra
pmde_g_vector <- hipparcos_gdr2$gaia_pmde
pmra_sd_g_vector <- hipparcos_gdr2$gaia_pmra_e
pmde_sd_g_vector <- hipparcos_gdr2$gaia_pmde_e
pmra_pmde_g_corr <- hipparcos_gdr2$gaia_pmra_pmde
pmra_var_g_vector <- pmra_sd_g_vector^2
pmde_var_g_vector <- pmde_sd_g_vector^2
pmra_pmde_g_cov <- pmra_pmde_g_corr*pmra_sd_g_vector*pmde_sd_g_vector

pmra_h_vector <- hipparcos_gdr2$pmRA
pmde_h_vector <- hipparcos_gdr2$pmDE
pmra_sd_h_vector <- hipparcos_gdr2$e_pmRA
pmde_sd_h_vector <- hipparcos_gdr2$e_pmDE
pmra_pmde_h_corr <- hipparcos_gdr2$pmDE.pmRA
pmra_var_h_vector <- pmra_sd_h_vector^2
pmde_var_h_vector <- pmde_sd_h_vector^2
pmra_pmde_h_cov <- pmra_pmde_h_corr*pmra_sd_h_vector*pmde_sd_h_vector


comp1_vector <- pmra_g_vector-pmra_h_vector
comp2_vector <- pmde_g_vector-pmde_h_vector

delta_p <- matrix(0, nrow = nrow(hipparcos_gdr2), ncol=2)


for (i in seq(1, length(comp1_vector))) {
  delta_p[i, 1]=comp1_vector[i]
  delta_p[i, 2]=comp2_vector[i]
}

delta_p <- cbind(comp1_vector, comp2_vector)
head(delta_p, 2)


C_g <- cbind(pmra_var_g_vector, pmde_var_g_vector, pmra_pmde_g_cov )
C_h <-  cbind(pmra_var_h_vector, pmde_var_h_vector, pmra_pmde_h_cov )
C <- C_g + C_h
head(C, 2)
colnames(C) <- c("pmra_v", "pmde_v", "pmra_pmde_cov")



muestra_p_m <- numeric(nrow(hipparcos_gdr2))


for (i in seq(1, nrow(hipparcos_gdr2))) {
  inv_det <- 1/(C[i,1]*C[i,2]-C[i, 3]^2)
  muestra_p_m[i]=delta_p[i,1]^2*inv_det*C[i,2]+delta_p[i, 1]*delta_p[i,2]*2*inv_det*(-C[i,3])+delta_p[i,2]^2*inv_det*C[i,1]
  
  
  
}

head(muestra_p_m)
length(muestra_p_m)

ks.test(muestra_p_m, 'pchisq', 2)

qqPlot(muestra_p_m, dist='chisq', df=2, ylim=c(0, 25), col='orange', xlab='Theoretical Quantiles', ylab='Sample Quantiles',
       main='Proper Motions', grid = F, envelope = 0.95)

df_muestra_pm <- as.data.frame(muestra_p_m)
length(df_muestra_pm$muestra_p_m)
ordered <- order(df_muestra_pm$muestra_p_m)
index <- ordered[seq(5000, 5701)]
auxiliary <- df_muestra_pm$muestra_p_m[-index]
df_auxiliary <- as.data.frame(auxiliary)
length(auxiliary)
ggplot(data.frame(df_auxiliary), aes(sample = df_auxiliary$auxiliary)) + 
  stat_qq(distribution = stats::qchisq, dparams = list(df = 2), col = 'orange') +
  stat_qq_line(distribution = stats::qchisq, dparams = list(df = 2), col = 'blue')

ranks <- rank(muestra_p_m, ties.method="max")
F_plot <- (ranks-0.5)/length(muestra_p_m)
F0inv <- qchisq(F_plot, df = 2)
df_plot <- data_frame(muestra_p_m, F0inv)
p <- ggplot(df_plot, aes(F0inv, muestra_p_m))
p + geom_point(col = 'blue') + geom_abline(intercept = 0, slope = 1, col = 'orange') + coord_cartesian(xlim =c(0, 5), ylim = c(0, 9)) +
  xlab('Theoretical quantiles') + ylab('Sample quantiles') + labs(title = 'Proper Motions')

# Graphics with ggplot 2

library(latex2exp)
a <- (omega_g-omega_h)^2
b <- abs(hipparcos_gdr2$deltaHp)
log_a <- log(a, 10)
model1 <- lm(b~log_a)
summary(model1)
cor(b, log_a)

library(ggplot2)
df <- data.frame(log_a, b)
plot <- ggplot(df, aes(b, log_a))
plot + geom_point()
plot + geom_point(alpha=1/2, size=2, color='blue') +
  geom_smooth(method='lm', col='green') + labs(x=TeX('$|\\Delta m_V|$')) + labs(y=TeX('$log(V)$')) +
  scale_y_continuous(limit = c(-3,3)) + ggtitle('Parallaxes') + theme_grey() + theme(legend.position = '')


almacenar <- numeric(nrow(hipparcos_gdr2))
for (i in seq(1, nrow(hipparcos_gdr2))) {
  almacenar[i]=delta_p[i, 1]^2+delta_p[i, 2]^2}

log_almacenar <- log(almacenar, 10)
plot(log_almacenar~b, ylab='log (W·W)', xlab='absolute value of visual magnitude difference', ylim=c(-3, 3), col='orange', main='Proper Motions')
model2 <- lm(log_almacenar~b)
abline(model2, col='blue')
summary(model2)
df1 <- data.frame(log_almacenar, b)
plot <- ggplot(df, aes(b, log_almacenar))
plot + geom_point()
plot + geom_point(alpha=1/2, size=2, color='blue') +
  geom_smooth(method='lm', col='green') + labs(x=TeX('$|\\Delta m_V|$')) + ggtitle('Proper Motions') +
  labs(y=TeX('$log(W·W)$')) + theme(legend.position ='')
cor(b, log_almacenar)

# Histograms

library(ggplot2)
library(gbm)
library(latex2exp)

head(hipparcos_gdr2)
max(hipparcos_gdr2$rho)
max(hipparcos_gdr2$Hp)
max(hipparcos_gdr2$Plx)

summary(hipparcos_gdr2$rho)
summary(hipparcos_gdr2$Hp)
summary(hipparcos_gdr2$Plx)

quantile(hipparcos_gdr2$rho, probs = c(0.05, 0.95))
quantile(hipparcos_gdr2$Hp, probs = c(0.05, 0.95))
quantile(hipparcos_gdr2$Plx, probs = c(0.05, 0.95))

par(mfrow = c(1,1))
hist1 <- ggplot(data = hipparcos_gdr2,
                aes(hipparcos_gdr2$rho)) + geom_histogram(col="red", fill="green")+
  labs(x=TeX("Angular separation (\\textit{arcsec})"), y="Number of double stars")

hist2 <- ggplot(data = hipparcos_gdr2,
                aes(hipparcos_gdr2$Hp)) + geom_histogram(col="red", fill="green")+
  labs(x=TeX("Visual magnitude, $m_V$ \\textit{(mag)}"), y="Number of double stars")

hist3 <- ggplot(data = hipparcos_gdr2,
                aes(hipparcos_gdr2$Plx)) + geom_histogram(col="red", fill="green", bins = 120)+
  labs(x=TeX("Parallax \\textit{(mas)}"), y="Number of double stars") + coord_cartesian(xlim = c(0, 40), ylim = c(0,750))

grid.arrange(hist1, hist2, hist3, nrow=3)

