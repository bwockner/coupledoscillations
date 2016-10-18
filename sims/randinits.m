%script to randomly vary the initial conditions of the sim
clear
theta0=deg2rad(70);
x0=0.15;
initconds=[x0,theta0];
ninits=10;
errs=[0.0005,deg2rad(1)]; %uncertainty of 0.5mm and 1 degree
inx=zeros(ninits);
intheta=zeros(ninits);
for i=1:ninits;
    %randomly add an error to the initial conditions
    rands=(rand(1,2)-0.5*ones(1,2));
    rands=rands.*errs;
    x=x0+rands(1);
    theta=theta0+rands(2);
    inx(i)=x;
    intheta(i)=theta;
    name=int2str(i);
    Project_iteration(x,theta,name);
end
dlmwrite('inits.txt', inx);
dlmwrite('inits.txt', intheta, '-append')