#The following code reproduce the FWER and FDR boxplots in Section 4.2.
set.seed(1)
p=0.01


par(mfrow=c(3, 3), mgp=c(2, 0.5, 0), mar=c(3, 5, 2, 1)+0.1)

for (t in 1:3) {
  p=0.01
  n=5*10^t
  theta=rbinom(n,1,p)
  x=(1-theta)*rnorm(n)+theta*rnorm(n,5,1)
  
  
  
  boxplot(x,main='Tukey-type Boxplot',outpch=19)
  
  mtext(text = paste0("k = ", t),
        side = 2,        # y-axis (left side)
        line = 3,        # distance from plot
        at = mean(range(par("usr")[3:4])),  # center vertically
        font = 2,        # bold
        cex = 1.1)
  # slightly larger text
  
  holm.boxplot(x,0.01)
  
  bh.boxplot(x,0.01)

  
}
