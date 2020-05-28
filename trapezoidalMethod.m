function [yNext] = trapezoidalMethod(t,Y)
  A = [0  0  1  0;
       0  0  0  1; 
       0  0  0  1 ; 
       0  0 -1  0 ];
  B = [0;0;1;0];
  I = eye(4);
  yNext = ((I-0.5*t*A)^(-1))*(Y+0.5*t*(A*Y+2*B));
end
