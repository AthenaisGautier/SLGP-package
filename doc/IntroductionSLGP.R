## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
collapse = TRUE,
comment = "#>"
)

## ----loadHousing, warning=FALSE-----------------------------------------------
library(dplyr)
# Load the dataset (available in MASS package)
if (!requireNamespace("MASS", quietly = TRUE)) install.packages("MASS")
data("Boston", package = "MASS")
df <- Boston %>%
  mutate(age_bin = cut(age, breaks = seq(0, 100, by = 10), include.lowest = FALSE)) %>%
  group_by(age_bin) %>%
  mutate(age_bin = paste0(age_bin, "\nn=", n()))%>%
  ungroup()%>%
  mutate(age_bin = factor(age_bin, 
                          levels = sort(unique(age_bin), decreasing = FALSE))) %>%
  data.frame()

range_response <- c(0, 50) # Can use range(df$medv), or user defined range as we do here
range_x <- c(0, 100) # Can use range(df$age), or user defined range as we do here

## ----figureHousing, fig.cap = "A visual representation of the dependency of the median value of owner-occupied homes on proportion of owner-occupied units constructed before 1940 in the Boston Housing dataset.", fig.fullwidth=TRUE, fig.height=4, fig.width=10, fig.align='center',fig.pos="H"----
library(ggplot2)
library(ggpubr)
library(viridis)
# Scatterplot: med vs. age
scatter_plot <- ggplot(df, aes(x = age, y = medv)) +
  geom_point(alpha = 0.5, color = "navy") +
  labs(x = "Proportion of owner-occupied units\nbuilt prior to 1940 [AGE, %]", 
       y = "Median value of owner-occupied homes [MEDV, k$]",
       title = "Median value vs. age of homes") +
  theme_bw()+
  coord_cartesian(xlim=range_x,
                  ylim=range_response)

# Histogram: Distribution of med by age bin
hist_plot <- ggplot(df, aes(x = medv)) +
  geom_histogram(mapping=aes(y=after_stat(density)),
                 position = "identity", breaks = seq(0, 50, 2.5),
                 fill="darkgrey", col="grey50", lwd=0.2, alpha=0.7) +
  geom_rug(sides = "b", color = "navy", alpha = 0.5)+
  facet_wrap(~ age_bin, scales = "free_y", nrow=2) +
  labs(x = "Median value of owner-occupied homes [MEDV, k$]", 
       y = "Probability density", 
       title = "Histogram of median housing values by AGE group") +
  theme_bw()+
  coord_cartesian(xlim=range_response,
                  ylim=c(0, 0.25))
ggarrange(scatter_plot, hist_plot, ncol = 2, nrow = 1,
          widths = c(0.3, 0.7))
# ggsave("./Figures/scatter.pdf", width=10, height=3.5)

## ----SLGPfitting00------------------------------------------------------------
library(SLGP)
modelMAP <- slgp(medv~age, # Use a formula to specify predictors VS response
                 # Can use medv~. for all variables,
                 # Or medv ~ age + var2 + var3 for more variables
                 data=df,
                 method="MAP", #Maximum a posteriori estimation scheme
                 basisFunctionsUsed = "RFF",
                 interpolateBasisFun="WNN", # Accelerate inference
                 hyperparams = list(lengthscale=c(0.15, 0.15), 
                                    # Applied to normalised data
                                    # So 0.15 is 15% of the range of values
                                    sigma2=1), 
                 # Will be re-selected with sigmaEstimationMethod
                 sigmaEstimationMethod = "heuristic", 
                 # Set to heuristic for numerical stability                 
                 predictorsLower= c(range_x[1]),
                 predictorsUpper= c(range_x[2]),
                 responseRange= range_response,
                 opts_BasisFun = list(nFreq=100,
                                      MatParam=5/2),
                 seed=1)

## ----SLGPplotting2, fig.cap = "Predictive probability density of medv at age, seen over slices.", fig.fullwidth=TRUE, fig.height=5, fig.width=10, fig.align='center',fig.pos="H"----


library(SLGP)
library(viridis)
dfGrid <- data.frame(expand.grid(seq(range_x[1], range_x[2], 5), 
                                 seq(range_response[1], range_response[2],, 101)))
colnames(dfGrid) <- c("age", "medv")
pred <- predictSLGP_newNode(SLGPmodel=modelMAP,
                            newNodes = dfGrid)

scale_factor <- 100
ggplot()  +
  labs(y = "Proportion of owner-occupied units\nbuilt prior to 1940 [%]", 
       x = "Median value of owner-occupied homes [k$]",
       title = "Median value vs. age of homes") +
  theme_bw()+
  geom_ribbon(data=pred,
              mapping=aes(x=medv, ymax=scale_factor*pdf_1+age, 
                          ymin=age, group=-age, fill=age),
              col="grey", alpha=0.9)+
  geom_point(data=df,
             mapping=aes(x = medv, y = age), alpha = 0.5, color = "navy")+
  scale_fill_viridis(option = "plasma",
                     guide = guide_colorbar(nrow = 1,
                                            title = "Indexing variable: Proportion of owner-occupied units built prior to 1940",
                                            barheight = unit(2, units = "mm"),
                                            barwidth = unit(55, units = "mm"),
                                            title.position = 'top',
                                            label.position = "bottom",
                                            title.hjust = 0.5))+
  theme(legend.position = "bottom")+
  coord_flip()


## ----SLGPfitting2-------------------------------------------------------------
modelLaplace <- retrainSLGP(SLGPmodel=modelMAP, 
                            newdata = df, 
                            method="Laplace")

## ----SLGPfitting22, eval=FALSE------------------------------------------------
#  # Or equivalent, more explicit in the re-using of the elements
#  # From the SLGP prior
#  
#  modelLaplace <- slgp(medv~age,
#                   data=df,
#                   method="MAP", #Maximum a posteriori estimation scheme
#                   basisFunctionsUsed = "RFF",
#                   interpolateBasisFun="WNN", # Accelerate inference
#                   hyperparams = modelMAP@hyperparams,
#                   sigmaEstimationMethod = "none",# Already selected in the prior
#                   predictorsLower= c(range_x[1]),
#                   predictorsUpper= c(range_x[2]),
#                   responseRange= range_response,
#                   opts_BasisFun = modelMAP@opts_BasisFun,
#                   BasisFunParam = modelMAP@BasisFunParam,
#                   seed=1)

## ----SLGPLaplaceplot, fig.cap = "Predictive probability density (and draws from a Laplace approximation) of medv at age, seen over 3 slices.", fig.fullwidth=TRUE, fig.height=4, fig.width=8, fig.align='center',fig.pos="H"----
# Define the three selected values
selected_values <- c(20, 50, 95)
gap <- 2.5

dfGrid <- data.frame(expand.grid(seq(3), 
                                 seq(range_response[1], range_response[2],, 101)))
colnames(dfGrid) <- c("ID", "medv")
dfGrid$age <- selected_values[dfGrid$ID]
pred <- predictSLGP_newNode(SLGPmodel=modelLaplace,
                            newNodes = dfGrid)
pred$meanpdf <- rowMeans(pred[, -c(1:3)])

library(tidyr)
# Filter the data: keep values within ±5 of the selected ones
df_filtered <- df %>%
  mutate(interval=findInterval(age, c(0, 
                                      selected_values[1]-gap, 
                                      selected_values[1]+gap, 
                                      selected_values[2]-gap, 
                                      selected_values[2]+gap, 
                                      selected_values[3]-gap, 
                                      selected_values[3]+gap)))%>%
  filter(interval %in% c(2, 4, 6))%>%
  group_by(interval)%>%
  mutate(category = paste0("Age close to ", c("", selected_values[1],
                                              "", selected_values[2],
                                              "", selected_values[3])[interval], 
                           "\nn=", n()))
names <- sort(unique(df_filtered$category))
pred$category <- names[pred$ID]

set.seed(1)
selected_cols <- sample(seq(1000), size=10, replace=FALSE)
df_plot <- pred %>%
  dplyr::select(c("age", "medv", "category", 
                  paste0("pdf_", selected_cols)))%>%
  pivot_longer(-c("age", "medv", "category"))


ggplot(mapping=aes(x = medv)) +
  geom_histogram(df_filtered,
                 mapping=aes(y=after_stat(density)),
                 position = "identity", breaks = seq(0, 50, 2.5),
                 fill="darkgrey", col="grey50", lwd=0.2, alpha=0.7) +
  geom_rug(data=df_filtered, sides = "b", color = "navy", alpha = 0.5)+
  geom_line(data=df_plot, mapping=aes(y=value, group=name), 
            color = "black", lwd=0.1, alpha=0.5)+
  geom_line(data=pred, mapping=aes(y=meanpdf, group=category), color = "red")+
  facet_wrap(~ category, scales = "free_y", nrow=1) +
  labs(x = "Median value of owner-occupied homes [k$]", 
       y = "Probability density", 
       title = "Histogram of median value at bins centered at several 'age' values, with width 5\nversus SLGP draws from a Laplace approximation (black curves) and average (red curve)") +
  theme_bw()+
  coord_cartesian(xlim=range_response,
                  ylim=c(0, 0.25))


## ----SLGPfitting3, eval=FALSE-------------------------------------------------
#  modelMCMC <- slgp(medv~age, # Use a formula to specify predictors VS response
#                    # Can use medv~. for all variables,
#                    # Or medv ~ age + var2 + var3 for more variables
#                    data=df,
#                    method="MCMC", #MCMC
#                    basisFunctionsUsed = "RFF",
#                    interpolateBasisFun="WNN", # Accelerate inference
#                    hyperparams = list(lengthscale=c(0.15, 0.15),
#                                       # Applied to normalised data
#                                       # So 0.15 is 15% of the range of values
#                                       sigma2=1),
#                    # Will be re-selected with sigmaEstimationMethod
#                    sigmaEstimationMethod = "heuristic", # Set to heuristic for numerical stability
#                    predictorsLower= c(range_x[1]),
#                    predictorsUpper= c(range_x[2]),
#                    responseRange= range_response,
#                    opts_BasisFun = list(nFreq=100,
#                                         MatParam=5/2),
#                    opts = list(stan_chains=2, stan_iter=1000))

## ----SLGPMCMCplotMoments, fig.cap = "Simultaneous prediction of the fields moments (and associated uncertainty) using a SLGP model", fig.fullwidth=TRUE, fig.height=4, fig.width=10, fig.align='center',fig.pos="H"----
dfX <- data.frame(age=seq(range_x[1], range_x[2], 2))

predMean <- predictSLGP_moments(SLGPmodel=modelLaplace,
                                newNodes = dfX, 
                                power=c(1),
                                centered=FALSE) # Uncentered moments
# For the mean
predVar <- predictSLGP_moments(SLGPmodel=modelLaplace,
                               newNodes = dfX, 
                               power=c(2), 
                               centered=TRUE) # Centered moments
# For the variance, Kurtosis and Skewness

pred <- rbind(predMean, predVar)
pred <- pred %>%
  pivot_longer(-c("age", "power"))%>%
  mutate(value=ifelse(power==2, sqrt(value), value))%>% # Define std
  data.frame()

pred$power <- factor(c("Expected value", 
                       "Standard deviation")[as.numeric(pred$power)],
                     levels=c("Expected value", "Standard deviation"))
df_plot <- pred %>%
  group_by(age, power)%>%
  summarise(q10 = quantile(value, probs=c(0.1)),
            q50 = quantile(value, probs=c(0.5)),
            q90 = quantile(value, probs=c(0.9)),
            mean = mean(value), .groups="keep")%>%
  ungroup() # summarise uncertainty


ggplot(df_plot, mapping=aes(x = age)) +
  geom_ribbon(mapping = aes(ymin=q10, ymax=q90),
              alpha = 0.25, lty=2, col="black", fill="cornflowerblue")+
  geom_line(mapping=aes(y=q50))+
  facet_grid(.~power)+
  labs(x = "Proportion of owner-occupied units built prior to 1940 [AGE, %]", 
       y = "Moment value") +
  theme_bw()+
  coord_cartesian(xlim=range_x,
                  ylim=range_response)


## ----SLGPMCMCplotQuantiles, fig.cap = "Simultaneous quantile prediction using a SLGP model, with estimation performed by MAP", fig.fullwidth=TRUE, fig.height=5, fig.width=10, fig.align='center',fig.pos="H"----
probsL <- c(5, 25, 50, 75, 95)/100
dfX <- data.frame(age=seq(range_x[1], range_x[2], 2))
pred <- predictSLGP_quantiles(SLGPmodel=modelLaplace,
                              newNodes = dfX, 
                              probs=probsL)

df_plot <- pred %>%
  pivot_longer(-c("age", "probs"))%>%
  group_by(age, probs)%>%
  summarise(q10 = quantile(value, probs=c(0.1)),
            q50 = quantile(value, probs=c(0.5)),
            q90 = quantile(value, probs=c(0.9)),
            mean = mean(value), .groups="keep")%>%
  ungroup()%>%
  mutate(probs=factor(paste0("Quantile: ", 100*probs, "%"),
                      levels=paste0("Quantile: ", 100*probsL, "%")))


plot1 <- ggplot(df_plot, mapping=aes(x = age)) +
  geom_ribbon(mapping = aes(ymin=q10, ymax=q90, 
                            col=probs, fill=probs, group=probs),
              alpha = 0.25, lty=3)+
  geom_line(mapping=aes(y=q50, lty="SLGP", 
                        col=probs,group=probs), lwd=0.75)+
  labs(x = "Proportion of owner-occupied units built prior to 1940 [AGE, %]", 
       y = "Median value of owner-occupied homes [MEDV, k$]",
       fill = "Quantile levels",
       col = "Quantile levels", 
       lty="Quantile estimation method") +
  theme_bw()+
  coord_cartesian(xlim=range_x,
                  ylim=range_response)+
  scale_linetype_manual(values=c(2, 1))

plot2 <- ggplot(df, mapping=aes(x=age))+
  geom_histogram(col="navy", fill="grey", bins = 30, alpha = 0.3,
                 breaks=seq(0, 100, 5))+
  theme_minimal()+
  labs(x = NULL, 
       y = "Sample count") 

library(ggpubr)
p_with_marginal <- ggarrange(
  plot2, plot1,
  ncol = 1, heights = c(1, 3),  # adjust height ratio
  align = "v"
)

print(p_with_marginal)

