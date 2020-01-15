load("timings.RData")
for (d in c(2,5,10))
{
  cat("\n")
  for (n in 4:6)
  {
    cat(mean(sapply(1:100, function(i) tm[[i]]$fm[d,n])))
    cat(" ")
  }
}
