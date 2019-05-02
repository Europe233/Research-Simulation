function [H] = H_matrix(index)
global fre n

time=(index-1)/fre;

h1=[sin(time);0;0];
h2=[0;3*cos(time);0];
h3=[0;0;2*sin(time)];
h4=[0;4*sin(time);0];

H=[h1;h2;h3;h4;zeros(2*n,1)];

end