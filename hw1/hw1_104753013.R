path = "pam1.txt"
#path = "/Users/blueswen/Documents/NCCU/1051/Bio/hw/hw1/pam1.txt"
# read PAM1 from data
pam1<-read.table(file = path)
pam1 <- as.matrix(pam1/10000)

d <- 250

result <- matrix()

for(i in 1:d){
  if(i == 1){
    result <- pam1
  }
  else{
    result <- result %*% pam1
  }
}

result <- data.frame(round(result %*% diag(1/colSums(result))*100))
colnames(result) <- rownames(result)

write.table(result,file='pam250.txt',sep='\t',quote=FALSE,col.names=NA)