function F = proj1E100_freqresp(G,C,L,b,w,Vg)

mag = [];

for iter = 1:length(w)
   
    omega = w(iter);
    
    A = G+1i*omega*C+(1/(1i*omega))*L ;
    x = A\b;
    
    mag = [mag, abs(x)];
end

V3 = mag(3,1:length(w));

F = 20*log10(V3/abs(Vg));