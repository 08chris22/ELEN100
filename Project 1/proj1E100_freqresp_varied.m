function [F, Qout] = proj1E100_freqresp_varied(QX,b,w,Vg,num)

F = []; Qout = [];

rng('default');

for iter = 1:num

    Q = proj1E100_variation(QX);
   
    Qout = [Qout; Q];
    
    G1 = [ ...
                (1)       (0)               (0); ...
                (-1/Q(1)) (1/Q(1) + 1/Q(2)) (-1/Q(2)); ...
                (0)       (-1/Q(2))         (1/Q(2))];

    G3 = [ 0,0,0;0,0,0;0,0,0 ];

    G2 = [ 0,0,0;0,Q(3),0;0,0,Q(4) ];
   
    FX = proj1E100_freqresp( G1,G2,G3,b,w,Vg);
   
    F = [F; FX];

end
