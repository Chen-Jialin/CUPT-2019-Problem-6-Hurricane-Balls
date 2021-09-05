clear;clc;

%T = [2.1,12.3,26.3,30.8];
%f = [30,20,10,0];

T = [5.1,9.6,13.9,19.7,27.3,34.9,43,48.2];
f = [35,30,25,20,15,10,5,0];

plot(T,f,'.','MarkerSize',20);
xlabel('时间t/s','Interpreter','tex');
ylabel('进动频率/Hz','Interpreter','tex');
grid on;