% Script written by:
% Zhuo Li (zhuol7@student.unimelb.edu.au)
% The University of Melbourne
 
clear all
close all
clc
 
%% RECORD THE VIDEO
VIDEO = 1;
 
%% SETUP INITIAL CONDITIONS
tspan = [0 10];
init = [30*pi/180; 5*pi/180; 0*pi/180; 0*pi/180; 5.5; 0; 1.0; 20.0];
options = odeset('RelTol',1e-7,'AbsTol',1e-7');
sol = ode45(@func ,tspan,init,options); 
 
%% EVAULATE THE SOLUTION
dt = 0.02;
t = tspan(1):dt:tspan(2);
X = deval(sol,t);
 
%% PLOT THE STATES
plot(t,[X(1,:);X(2,:);X(3,:);X(4,:)],'LineWidth', 1.5)
xlabel('time', 'FontSize',14)
ylabel('states', 'FontSize',14)
h_1 = legend('$\alpha$','$\beta$','$\gamma$','$\delta$');
set(h_1,'Interpreter','latex')
savefig('fig1')
plot(t,[X(5,:);X(6,:);X(7,:);X(8,:)],'LineWidth', 1.5)
xlabel('time', 'FontSize',14)
ylabel('states', 'FontSize',14)
h_2 = legend('$\dot{\alpha}$','$\dot{\beta}$','$\dot{\gamma}$','$\dot{\delta}$');
set(h_2,'Interpreter','latex')
savefig('fig2')

%% SETUP PARAMETERS
r = 1;    	% Minor radius of torus frames and radius of cylindrical frame
R_f = 34.50;  % Major radius of torus frames
R_r = 28.50; 	% Radius of rotor
H_f = 94;   	% Height of cylindrical frame
H_r = 2;      	% Height of rotor
L = 47;  	% Distance from the tip to the mid-length of cylindrical frame
 
%% CREATE GYROSCOPE

% Create the rotor

theta = [0:10:360];

[x1,y1,z1] = cylinder(1);
z1 = z1*H_r;
x1 = x1*R_r;
y1 = y1*R_r;
z1 = z1+L;

x1s = R_r.*cosd(theta);
y1s = R_r.*sind(theta);
z1s = ones(1,length(theta)) .* z1(1);
 
x2s = x1s;
y2s = y1s;
z2s = z1s + H_r;
 
xp = [0, x2s(1)];
yp = [0, y2s(1)];
zp = [L + H_r, z2s(1)];
 
% Create the Frames
 
% Horizontal ring
xh = R_f .* cosd(theta);
yh = R_f .* sind(theta);
zh = L*ones(length(theta))+r;
 
% Vertical ring
xv = R_f.* cosd(theta);
yv = 0.*ones(length(theta));
zv = R_f .* sind(theta)+L+r;
 
% Cylinder
[xa,ya,za] = cylinder(0.5,30);
xa = xa;
ya = ya;
za = za*H_f;
 
%% SETUP VIDEO IF REQUIRED
if VIDEO
	FPS = 1/dt;
	MyVideo = VideoWriter('Gyroscope_simulation','MPEG-4');
	MyVideo.FrameRate = FPS;
	open(MyVideo);
end
 
%% CREATE ANIMATION
handle = figure;
hold on
for i = 1:length(t)
	cla
	a_0 = X(1,i);
	b_0 = X(2,i);
	c_0 = X(3,i);
	d_0 = X(4,i);
	a_1 = X(5,i);
	b_1 = X(6,i);
	c_1 = X(7,i);
	d_1 = X(8,i);
	
	% Redefine the rotation matrices required for simulation
	R01 = [cos(a_0) -sin(a_0) 0; sin(a_0) cos(a_0) 0; 0 0 1];
    R12 = [1 0 0; 0 cos(b_0) -sin(b_0); 0 sin(b_0) cos(b_0)];
    R23 = [cos(c_0) -sin(c_0) 0; sin(c_0) cos(c_0) 0; 0 0 1];
    R34 = [cos(d_0) -sin(d_0) 0; sin(d_0) cos(d_0) 0; 0 0 1];
	R03 = R01 * R12 * R23;
	R04 = R03 * R34;
	
	% Rotation frame {4} is defined so that it is attached to the rotor
	[x1_rotated, y1_rotated, z1_rotated] = rotation(x1,y1,z1,R04);
	[x1s_rotated, y1s_rotated, z1s_rotated] = rotation(x1s,y1s,z1s,R04);
	[x2s_rotated, y2s_rotated, z2s_rotated] = rotation(x2s,y2s,z2s,R04);
	[xp_rotated, yp_rotated, zp_rotated] = rotation(xp,yp,zp,R04);
	
	% Rotation frame {3} is attached to the gyroscope frame
	[xa_rotated, ya_rotated, za_rotated] = rotation(xa,ya,za,R03);
	[xv_rotated, yv_rotated, zv_rotated] = rotation(xv,yv,zv,R03);
	[xh_rotated, yh_rotated, zh_rotated] = rotation(xh,yh,zh,R03);
	
	% Plots the points / shapes defined
	surf(x1_rotated, y1_rotated, z1_rotated, 'EdgeColor','none');
	surf(xa_rotated, ya_rotated, za_rotated);
	plot3(xv_rotated, yv_rotated, zv_rotated, 'k', 'LineWidth', 3);
	plot3(xh_rotated, yh_rotated, zh_rotated, 'k', 'LineWidth', 3);
	fill3(x1s_rotated, y1s_rotated, z1s_rotated, 'g');
	fill3(x2s_rotated, y2s_rotated, z2s_rotated, 'g');
	plot3(xp_rotated, yp_rotated, zp_rotated, 'r', 'LineWidth', 2);
 
	axis square
	view(3)
	axis(1*[-50 50 -50 50 0 100])
	xlabel('X')
	ylabel('Y')
	zlabel('Z')
	if VIDEO
    	writeVideo(MyVideo,getframe(handle));
	else
    	pause(dt)
	end
end
 
if VIDEO
close(MyVideo)
end
