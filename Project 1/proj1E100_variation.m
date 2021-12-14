function Q = proj1E100_variation(R)

m = length(R) ;
g = ones(1,m)+0.5*rand(1,m) ;

Q = 0.8*(R.*g) ;