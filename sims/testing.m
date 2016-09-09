%Initial scratchpad for testing and things. probably not useful
%for anything else
clear
%%%Ball variables (might turn this into a class for convenience later)
Xball=zeros(3,2); %2*3-vector representing postion of ball COM in R3 and ball rotation
mball=1; %kg
rball=1; %radius of the ball (m)
Iball=eye(3)*mball*rball^2*2/5; %Inertia matrix of the ball (for simplicity assume that the 
%ball has no product (off diagonal) terms; inertia of a sphere, kg m^2
%Ixx=Iyy=Izz=2/5*m*r^2
Vball=zeros([3,2]);
%Total energy of the ball (kj)
Eball=0.5*(transpose(Vball(:,2))*Iball*Vball(:,2))...%kinetic energy due to spin
            +0.5*mball*transpose(Vball(:,1))*Vball(:,1); %translational kinetic energy
