function [dist] = getOffset(initialSpeed)
y0 = [1;0;initialSpeed;0];
options = odeset('Events', @event);
intervall = [-2 2];
y0 = [1;0;initialSpeed;0];
[t,y] = ode45(@func, intervall, y0, options);
exitPoint = [y(end,1); y(end,2)];
dist = exitPoint - [0;1];
end