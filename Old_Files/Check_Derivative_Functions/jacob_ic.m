function [h] = jacob_ic(a)
h(1,1)=0;
h(1,2)=a(6);
h(1,3)=0;
h(1,4)=0;
h(1,5)=a(5);
