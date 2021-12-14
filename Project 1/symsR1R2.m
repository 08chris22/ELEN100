syms R1 R2
c1 = .1e-6;
c2 = .1e-6;
a = 3000;
b = 20000;
eq1 = 1/(a*b) == R1*R2*c1*c2;
eq2 = (a+b)/(a*b) == R1*c1 + R2*c2 + R1*c2;
sol = solve([eq1,eq2],[R1,R2]);
R1 = double(sol.R1);
R2 = double(sol.R2)