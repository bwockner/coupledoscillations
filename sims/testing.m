%Initial scratchpad for testing and things. probably not useful
%for anything else
clear
%global things:
gravity=[0,0,-9.81]; %direction and magnitude of gravitational acceleration
%%%Ball variables (might turn this into a class for convenience later)
xball=zeros(3,2); %2*3-vector representing postion of ball COM in R3 and ball rotation
mball=1; %kg
rball=1; %radius of the ball (m)
Iball=eye(3)*mball*rball^2*2/5; %Inertia matrix of the ball (for simplicity assume that the 
%ball has no product (off diagonal) terms; inertia of a sphere, kg m^2
%Ixx=Iyy=Izz=2/5*m*r^2
vball=zeros([3,2]);
%Total energy of the ball (kj)
Tball=0.5*(transpose(vball(:,2))*Iball*vball(:,2))...%kinetic energy due to spin
            +0.5*mball*transpose(vball(:,1))*vball(:,1); %translational kinetic energy

Fball=zeros([3,2]); %forces and moments acting on the ball
Fball(:,1)=mball*gravity; 

Vball= -mball*dot(gravity,xball(:,1));
%Velocity of the ball is tangential to the track:
%generally, this means that the ball's velocity must be tangent to the
%vector space that formed by the solution to u x v= 0. The point of contact
%between the track and the ball has to lie in this plane too
%
%Track:
%The system will act as a pendulum, where the centre of mass of the track
%swings below its centre of rotation. Lagrange's equations can be used to
%determine the equations oif motion for the ball and track.