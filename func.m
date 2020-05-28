function [yPrim] = func(t,Y)
A = [0  0  1  0;
     0  0  0  1; 
     0  0  0  1 ; 
     0  0 -1  0 ];
B = [0;0;1;0];
yPrim = A*Y + B;
% y' =[u;v;1+v;-u]
end