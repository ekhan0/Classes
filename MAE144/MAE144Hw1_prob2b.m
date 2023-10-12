%Define all the poles of f where we want them.
f_poles=[-1,-1,-3,-3,-6,-6];

%Define poles for a and zeros for b based on G(s) we are given. 
a_poles = [-1 1 -3 3 -6 6];
b_zeros = [-2 2 -5 5];

%Define all our polynomials based on our known poles and zeros.
f = RR_poly(f_poles,1);
a = RR_poly(a_poles,1);
b = RR_poly(b_zeros,1);

%Get initial x and y from RR_diophantine
[x,y] = RR_diophantine(a,b,f);

%Count keeps track of how many poles at s = -20 we're adding.
count = 0;

%Loops till order of y is less than order of x. 
while x.n <= y.n
    f_poles(end+1) = -20;
    f = RR_poly(f_poles,1);
    count = count+1;
    [x,y] = RR_diophantine(a,b,f);
end

%Display new x and y.
x
y

%Make our D transfer function based on the new y and x we have. 
D = RR_tf(y,x)

%Display how many K's were added 
count

%Test to make sure that a*x+b*y = f is actually true. Want residual = 0.
test=trim(a*x+b*y);
residual=norm(f-test)