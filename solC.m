%% (1) Euler forward
clear; clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% AB
unitCircle(0,0,1);
xlim([-1.5 1.5])
ylim([-1.5 1.5])
hold on

x = 1; u = -2;  
y = 0; v = 0;
tspan = [-2 2];
[t,Y] = ode45(@func, tspan, [1 0 -2 0]);

% plot(Y(:,1), Y(:,2)), old plot

% graph points
Xpoints = Y(:,1);
Ypoints = Y(:,2);

% circle points
r = 1; x = 0; y = 0;
th = (pi/2):pi/100:pi;
xunit = r * cos(th) + x;
yunit = r * sin(th) + y;

% amount of points on plot
sizeXvector = size(Ypoints(:,1));
pointsFunc = sizeXvector(1,1);

sizeXunitVector = size(xunit);
pointsOnCircle = sizeXunitVector(1,2);

tol = 0.013;
toInvesValX = [];
toInvesValY = [];
for i = 1:pointsFunc
    for j = 1:pointsOnCircle
        if (abs(Xpoints(i,1) - xunit(1,j)) < tol) && (abs(Ypoints(i,1) - yunit(1,j)) < tol)
            toInvesValX = [toInvesValX, Xpoints(i,1)];
            toInvesValY = [toInvesValY, Ypoints(i,1)]; 
        end
    end
end

interSec = [toInvesValX(1,1), toInvesValY(1,1)];
indexX = find(Xpoints == interSec(1,1));
indexY = find(Ypoints == interSec(1,2));

newY = [];
newX = [];
if (indexX == indexY)
    pointsInNewPlot = indexX;
    for i = 1:pointsInNewPlot
    newY = [newY, Ypoints(i,1)];
    newX = [newX, Xpoints(i,1)];
    end
end

hold on 
plot(newX, newY)

% choose speed
uPrelim = Y(:,3);
vPrelim = Y(:,4);
u = uPrelim(indexX,1);
v = vPrelim(indexY,1); % indexX=indexY=(33 in AB)

sizeSpeedVector = size(uPrelim);
sizeSpeed = sizeSpeedVector(1,1);

% fill x and y vectors with coordinates for straight line using u(33,1),
% v(33,1)
t = 0.1;
i = pointsInNewPlot;
straightX = [];
straightY = [];
x = Xpoints(indexX,1);
y = Ypoints(indexY,1);
while i ~= (sizeSpeed)
   straightX =[straightX ,x];
   straightY =[straightY ,y];
   x = x + u*t;
   y = y + v*t;
   i = i + 1;
end

% x = x + u33*h
% y = y + v33*h
for j = 1:(i-pointsInNewPlot-1)
   hold on
   plot([straightX(1,j) straightX(1,j+1)], [straightY(1,j) straightY(1,j+1)], 'Color', [0.3010, 0.7450, 0.9330]); 
end

for j = 1:pointsInNewPlot-1
   hold on
   plot([Xpoints(j,1) Xpoints(j+1,1)], [Ypoints(j,1) Ypoints(j+1,1)], 'Color', [0.3010, 0.7450, 0.9330]); 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% AB ends, 1) Euler forwards begins

unitCircle(0,0,1);
xlim([-1.5 1.5])
ylim([-1.5 1.5])
hold on

n = 60; i = 0;
Y = [1 0 -2 0];
t = 0.1;
xCoor = [];
yCoor = [];
xSpeed = [];
ySpeed = [];
while i~=n
   xCoor = [xCoor, Y(1,1)];
   yCoor = [yCoor, Y(1,2)];
   xSpeed = [xSpeed, Y(1,3)];
   ySpeed = [ySpeed, Y(1,4)];
   nextY = eulerForward(t,Y);
   Y = nextY; 
   i = i + 1; 
end

% % first plot
% n = 60; i = 1; 
% while i~=n
% plot([xCoor(1,i) xCoor(1,i+1)],[yCoor(1,i) yCoor(1,i+1)], 'Color', 'r');    
% hold on    
% i = i + 1;    
% end

% circle points
r = 1; x = 0; y = 0;
th = (pi/2):pi/200:pi;
xunit = r * cos(th) + x;
yunit = r * sin(th) + y;

% amount of points on plot
sizeXvector = size(yCoor(1,:));
pointsFunc = sizeXvector(1,2);

sizeXunitVector = size(xunit);
pointsOnCircle = sizeXunitVector(1,2);

%find points of intersection
tol = 0.13;
toInvesValX = [];
toInvesValY = [];
for i = 1:pointsFunc
    for j = 1:pointsOnCircle
        if (abs(xCoor(1,i) - xunit(1,j)) < tol) && (abs(yCoor(1,i) - yunit(1,j)) < tol)
            toInvesValX = [toInvesValX, xCoor(1,i)];
            toInvesValY = [toInvesValY, yCoor(1,i)]; 
        end
    end
end

interSec = [toInvesValX(1,1), toInvesValY(1,1)];
indexX = find(xCoor == interSec(1,1));
indexY = find(yCoor == interSec(1,2));

% save points up and until intersection, use speed of last point, make line
newY = [];
newX = [];
if (indexX == indexY)
    pointsInNewPlot = indexX;
    for i = 1:pointsInNewPlot
    newY = [newY, yCoor(1,i)];
    newX = [newX, xCoor(1,i)];
    end
end

% % choose speed
u = xSpeed(1,indexX);
v = ySpeed(1,indexY); % indexX=indexY=(33 in AB)

sizeSpeedVector = size(xCoor);
sizeSpeed = sizeSpeedVector(1,2);

% fill x and y vectors with coordinates for straight line using u(33,1),
% v(33,1)
t = 0.1;
i = pointsInNewPlot;
straightX = [];
straightY = [];
x = xCoor(1,indexX);
y = yCoor(1,indexY);
while i ~= (sizeSpeed)
   straightX =[straightX ,x];
   straightY =[straightY ,y];
   x = x + u*t;
   y = y + v*t;
   i = i + 1;
end

% x = x + u33*h
% y = y + v33*h
for j = 1:(i-pointsInNewPlot-1)
   hold on
   plot([straightX(1,j) straightX(1,j+1)], [straightY(1,j) straightY(1,j+1)], 'Color', [0.6350, 0.0780, 0.1840]); 
end

n = indexX; i = 1; 
if indexX == indexY
    while i~=n
    plot([xCoor(1,i) xCoor(1,i+1)],[yCoor(1,i) yCoor(1,i+1)], 'Color', [0.6350, 0.0780, 0.1840]);    
    hold on    
    i = i + 1;    
    end
end

% %% (2) Euler Backward
% clear; clc;
% 
% unitCircle(0,0,1);
% xlim([-1.5 1.5])
% ylim([-1.5 1.5])
% hold on

n = 60; i = 0;
Y = [1; 0; -2; 0];
t = 0.1;
xCoor = [];
yCoor = [];
xSpeed = [];
ySpeed = [];
while i~=n
   xCoor = [xCoor, Y(1,1)];
   yCoor = [yCoor, Y(2,1)];
   xSpeed = [xSpeed, Y(3,1)];
   ySpeed = [ySpeed, Y(4,1)];
   nextY = eulerBackward(t,Y);
   Y = nextY; 
   i = i + 1; 
end

% first plot
% n = 60; i = 1;
% while i~=n
% plot([xCoor(1,i) xCoor(1,i+1)],[yCoor(1,i) yCoor(1,i+1)], 'Color', 'g');    
% hold on    
% i = i + 1;    
% end

% circle points
r = 1; x = 0; y = 0;
th = (pi/2):pi/200:pi;
xunit = r * cos(th) + x;
yunit = r * sin(th) + y;

% amount of points on plot
sizeXvector = size(yCoor(1,:));
pointsFunc = sizeXvector(1,2);

sizeXunitVector = size(xunit);
pointsOnCircle = sizeXunitVector(1,2);

%find points of intersection
tol = 0.01;
toInvesValX = [];
toInvesValY = [];
for i = 1:pointsFunc
    for j = 1:pointsOnCircle
        if (abs(xCoor(1,i) - xunit(1,j)) < tol) && (abs(yCoor(1,i) - yunit(1,j)) < tol)
            toInvesValX = [toInvesValX, xCoor(1,i)];
            toInvesValY = [toInvesValY, yCoor(1,i)]; 
        end
    end
end

interSec = [toInvesValX(1,1), toInvesValY(1,1)];
indexX = find(xCoor == interSec(1,1));
indexY = find(yCoor == interSec(1,2));

% save points up and until intersection, use speed of last point, make line
newY = [];
newX = [];
if (indexX == indexY)
    pointsInNewPlot = indexX;
    for i = 1:pointsInNewPlot
    newY = [newY, yCoor(1,i)];
    newX = [newX, xCoor(1,i)];
    end
end

% % choose speed
u = xSpeed(1,indexX);
v = ySpeed(1,indexY); % indexX=indexY=(33 in AB)

sizeSpeedVector = size(xCoor);
sizeSpeed = sizeSpeedVector(1,2);

% fill x and y vectors with coordinates for straight line using u(33,1),
% v(33,1)
t = 0.1;
i = pointsInNewPlot;
straightX = [];
straightY = [];
x = xCoor(1,indexX);
y = yCoor(1,indexY);
while i ~= (sizeSpeed)
   straightX =[straightX ,x];
   straightY =[straightY ,y];
   x = x + u*t;
   y = y + v*t;
   i = i + 1;
end

% x = x + u33*h
% y = y + v33*h
for j = 1:(i-pointsInNewPlot-1)
   hold on
   plot([straightX(1,j) straightX(1,j+1)], [straightY(1,j) straightY(1,j+1)], 'Color', [0.4940, 0.1840, 0.5560]); 
end

n = indexX; i = 1; 
if indexX == indexY
    while i~=n
    plot([xCoor(1,i) xCoor(1,i+1)],[yCoor(1,i) yCoor(1,i+1)], 'Color', [0.4940, 0.1840, 0.5560]);    
    hold on    
    i = i + 1;    
    end
end
% %% (3) Trapezoidal method
% clear; clc;
% 
% unitCircle(0,0,1);
% xlim([-1.5 1.5])
% ylim([-1.5 1.5])
% hold on

n = 60; i = 0;
Y = [1; 0; -2; 0];
t = 0.1;
xCoor = [];
yCoor = [];
xSpeed = [];
ySpeed = [];
while i~=n
   xCoor = [xCoor, Y(1,1)];
   yCoor = [yCoor, Y(2,1)];
   xSpeed = [xSpeed, Y(3,1)];
   ySpeed = [ySpeed, Y(4,1)];
   nextY = trapezoidalMethod(t,Y);
   Y = nextY; 
   i = i + 1; 
end

% first plot
% n = 60; i = 1;
% while i~=n
% plot([xCoor(1,i) xCoor(1,i+1)],[yCoor(1,i) yCoor(1,i+1)], 'Color', 'b');    
% hold on    
% i = i + 1;    
% end

% circle points
r = 1; x = 0; y = 0;
th = (pi/2):pi/200:pi;
xunit = r * cos(th) + x;
yunit = r * sin(th) + y;

% amount of points on plot
sizeXvector = size(yCoor(1,:));
pointsFunc = sizeXvector(1,2);

sizeXunitVector = size(xunit);
pointsOnCircle = sizeXunitVector(1,2);

%find points of intersection
tol = 0.03;
toInvesValX = [];
toInvesValY = [];
for i = 1:pointsFunc
    for j = 1:pointsOnCircle
        if (abs(xCoor(1,i) - xunit(1,j)) < tol) && (abs(yCoor(1,i) - yunit(1,j)) < tol)
            toInvesValX = [toInvesValX, xCoor(1,i)];
            toInvesValY = [toInvesValY, yCoor(1,i)]; 
        end
    end
end

interSec = [toInvesValX(1,1), toInvesValY(1,1)];
indexX = find(xCoor == interSec(1,1));
indexY = find(yCoor == interSec(1,2));

% save points up and until intersection, use speed of last point, make line
newY = [];
newX = [];
if (indexX == indexY)
    pointsInNewPlot = indexX;
    for i = 1:pointsInNewPlot
    newY = [newY, yCoor(1,i)];
    newX = [newX, xCoor(1,i)];
    end
end

% % choose speed
u = xSpeed(1,indexX);
v = ySpeed(1,indexY); % indexX=indexY=(33 in AB)

sizeSpeedVector = size(xCoor);
sizeSpeed = sizeSpeedVector(1,2);

% fill x and y vectors with coordinates for straight line using u(33,1),
% v(33,1)
t = 0.1;
i = pointsInNewPlot;
straightX = [];
straightY = [];
x = xCoor(1,indexX);
y = yCoor(1,indexY);
while i ~= (sizeSpeed)
   straightX =[straightX ,x];
   straightY =[straightY ,y];
   x = x + u*t;
   y = y + v*t;
   i = i + 1;
end

% x = x + u33*h
% y = y + v33*h
for j = 1:(i-pointsInNewPlot-1)
   hold on
   plot([straightX(1,j) straightX(1,j+1)], [straightY(1,j) straightY(1,j+1)], 'Color', [0.8500, 0.3250, 0.0980]); 
end

n = indexX; i = 1; 
if indexX == indexY
    while i~=n
    plot([xCoor(1,i) xCoor(1,i+1)],[yCoor(1,i) yCoor(1,i+1)], 'Color', [0.8500, 0.3250, 0.0980]);    
    hold on    
    i = i + 1;    
    end
end