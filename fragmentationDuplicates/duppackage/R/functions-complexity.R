
libCompNewParam <-
  function(X,R,N){
    (log(N*(R-1)+X)-log(X)+N/X)^2
  }
