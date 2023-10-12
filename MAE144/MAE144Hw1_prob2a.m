%Could have written this in 2 lines like said on the homework,
%but I prefer the well broken down code so it is easy for others to follow.

%Define all the poles of f where we want them.
f_poles=[-1,-1,-3,-3,-6,-6];

%Define poles for a and zeros for b based on G(s) we are given. 
a_poles = [-1 1 -3 3 -6 6];
b_zeros = [-2 2 -5 5];

%Define all our polynomials based on our known poles and zeros.
f = RR_poly(f_poles,1);
a = RR_poly(a_poles,1);
b = RR_poly(b_zeros,1);

%Get x and y from RR_diophantine
[x,y] = RR_diophantine(a,b,f)

%Make our D transfer function based on the y and x we have. 
D = RR_tf(y,x)

%Test to make sure that a*x+b*y = f is actually true. Want residual = 0.
test=trim(a*x+b*y); 
residual=norm(f-test)