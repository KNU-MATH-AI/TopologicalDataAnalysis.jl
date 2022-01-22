X = [1 2 2 4; 1 2 6 1; 7 8 9 0]
rswap!(X,1,2)
smith(X)
gauss!(X)

X = X[2:end,2:end]