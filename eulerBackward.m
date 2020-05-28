function [yNext] = eulerBackward(t,Y)
  A = [0  0  1  0;
       0  0  0  1; 
       0  0  0  1 ; 
       0  0 -1  0 ];
  B = [0;0;1;0];
  I = eye(4);
  yNext = ((I-t*A)^(-1))*(Y+t*B);
end