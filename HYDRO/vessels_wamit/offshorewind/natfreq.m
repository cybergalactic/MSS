function F = natfreq(x,m,k,w)
F = x - sqrt(k/interp1(w,m,x));