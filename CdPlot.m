x1=[10	20	40	47	50	80	100	150	300];
y1=[2.81	2.04	1.54	1.46	1.46	1.4	1.39	1.37	1.27];
x2=[20 80 150 300];
y2= [2.4 1.54 1.44 1.38];
x3 = [20 40 80 300];
y3 = [2.03	1.52 1.37 1.38];
x4=[10	20	40	80	100	300];
yp= [2.78	2.01	1.51	1.35	1.33	1.37];
x5=[10	20	40	50	80	100	150	300];
ys =[2.67	2.08 1.73 1.65	1.51 1.45 1.36 1.22];
x6=[20 40 80];
yt = [2.22 1.48 1.29];
figure
plot(x2,y2,'K-o',x1,y1,'g-*',x3,y3,'r-*',x4,yp,'b-*',x5,ys,'c-*',x6,yt,'m-*');
legend("Present", "Paper", "Ye et al", 'Park', "S&B", "Tretton")
xlabel('Re');
ylabel('Cd');