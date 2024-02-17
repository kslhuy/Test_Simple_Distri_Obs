function [x_dot,y1,y2,y3,y4] = system_fcn(A,H1,H2,H3,H4,x)



x_dot = A*x;


y1 = H1*x;
y2 = H2*x;
y3 = H3*x;
y4 = H4*x;
end

