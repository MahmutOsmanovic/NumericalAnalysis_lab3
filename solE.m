clear; clc;
% from the textbook, page 163, we find:
% y'= Ay + b is stable for small h iff y' = cy stable for h and every
% eigenvalue c to A.
  A = [0  0  1  0;
       0  0  0  1; 
       0  0  0  1; 
       0  0 -1  0];
  E = eig(A);
  h = 0.1;
  % Eulers forward method is stable if there for every eigenvalue c to A we
  % have that: abs(1+hc)<=1, implicit method (needs less computation)
  % Is not "a stabil"
  disp('Euler forward')
  disp(prod(abs(1+h*E)<=1))
  
  % Euler backwards is stable if it for every eigenvalue c to A we have
  % that abs(1-hc)>=1, explicit method (needs more computation)
  % "Is "A stabil"
  disp('Euler backward')
  disp(prod(abs(1-h*E)>=1))
  
  % The trapezoidal method is stable if it for every eigenvalue c to A we have
  % that abs(1+hc/2)/abs(1-hc/2)<=1, fexplicit method (needs more computation)
  % Is "A stabil"
  disp('Trapezoidal method')
  disp(prod(abs(1+(h/2)*E)./abs(1-(h/2)*E)<=1))
  
  disp(E)
  % stable if all Re(eigenvalues) < 0. Has one eigenvalue >0, not stable
  % "styva problem", use implicit methods, enables larger stepsize