function [Ynext] = eulerForward(t,Y)
  A = [0  0  1  0;
       0  0  0  1; 
       0  0  0  1 ; 
       0  0 -1  0 ];
  B = [0;0;1;0];
  transY = transpose(Y);
  yPrim = A*transY + B;
  Ynext = Y + transpose(yPrim)*t;
end

% y(t+h) = y(t)   + y'(t)h <=>
% y(t+h) = y(t) + f(t,y(t))h
% y_k+1 = y_k + f(t_k,y_k)h_k