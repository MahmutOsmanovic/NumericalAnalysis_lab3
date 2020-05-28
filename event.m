function[value,stop,direction]=event(t,y)
y1 = y(1);
y2 = y(2);
value=norm([y1(end),y2(end)]) - 1;
stop=true; 
direction=1; 
end
