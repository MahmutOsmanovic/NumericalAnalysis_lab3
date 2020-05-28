clear; clc;

t = linspace(0,1);
plot(cos(2*pi*t), sin(2*pi*t), 'Color', [0, 0.4, 0.8]);
xlim([-1.5 1.5])
ylim([-1.5 1.5])
grid on
hold on

x = 1; u = -2;  
y = 0; v = 0;
tspan = [-2 2];
[t,Y] = ode45(@func, tspan, [x y u v]);

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
   plot([straightX(1,j) straightX(1,j+1)], [straightY(1,j) straightY(1,j+1)], 'Color', [0.9, 0, 0.2]); 
end

for j = 1:pointsInNewPlot-1
   hold on
   plot([Xpoints(j,1) Xpoints(j+1,1)], [Ypoints(j,1) Ypoints(j+1,1)], 'Color', [0.9, 0, 0.2]); 
end