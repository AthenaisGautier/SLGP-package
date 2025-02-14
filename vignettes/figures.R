library(ggplot2)
library(viridis)

# Compute Legendre polynomials up to degree n for a vector of x values
legendre_table <- function(n, x_vals) {
  x_vals <- as.vector(x_vals)  # Ensure x_vals is a vector
  m <- length(x_vals)          # Number of x values
  P <- matrix(0, nrow = m, ncol = n + 1)  # Matrix to store results
  
  # Base cases
  P[, 1] <- 1                # P_0(x) = 1
  if (n > 0) P[, 2] <- x_vals # P_1(x) = x
  
  # Compute P_n(x) using recurrence relation
  for (k in 2:n) {
    P[, k + 1] <- ((2 * k - 1) * x_vals * P[, k] - (k - 1) * P[, k - 1]) / k
  }
  
  colnames(P) <- paste0("P_", 0:n)  # Label columns for readability
  return(P)
}


# Example: Evaluate Bernstein basis polynomials of degree n = 5 at x = 0.3
x <- seq(-1, 1,, 20001)
M <- matrix(nrow=3*3, ncol=length(x))
for(i in seq(3)){
  p <- c(10, 100, 1000)[i]
  B <- legendre_table(p-1, x)
  for(j in seq(3)){
    set.seed(j)
    epsilon <- rnorm(p)
    
    f <-   B%*%epsilon
    SLGP <- f-max(f)
    SLGP <- SLGP *5 / diff(range(SLGP))
    SLGP <- exp(SLGP)/mean(exp(SLGP))
    M[(i-1)*3+j, ] <- SLGP
  }
} 
p <- c(sapply(c(10, 100, 1000), function(i){rep(i, (3*length(x)))}))
d <- c(sapply(seq(3), function(i){rep(i, length(x))}))
df <- data.frame(x=x, y=c(t(M)),
                 p=paste0("Number of basis functions:", p),
                 draw=paste0("Draw from the prior n°", d))
ggplot(df) +
  geom_line(mapping=aes(x=x, y=y, col=draw, group=draw))+
  facet_wrap(.~p, scale="free_y") +
  theme_bw()+
  theme(legend.position="bottom",
        legend.direction = "horizontal", 
        legend.title = element_blank())+
  coord_cartesian(ylim=c(0, 5))+
  xlab("Response")+
  ylab("Probability density")

ggsave(paste0("./Figures/badBasisFun",  ".pdf"), width=10, height=4.5)

