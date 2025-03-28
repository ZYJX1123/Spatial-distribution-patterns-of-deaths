close all;
clear; clc;

n=input( 'please input number of points n=');
pic=imread('');%Type in the file name of the image

I=imshow(pic);

loc_points=zeros(n,2);

for i=1:1:n
    
hold on;  
[x, y]=ginput(1);

hold on;
plot(x,y,'r.')
 
loc_points(i,1) = x;
loc_points(i,2) = y;

str=['  X:' num2str(x') ', Y:' num2str(y')];

end
