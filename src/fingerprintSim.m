% fingerprintSim.m 
% Kevin Coltin 
% 
% Solves a reaction diffusion equation for fingerprint formation. 
% Solves an equation of the form 
%   dc/dt = f(c) + D del^2 c, 
% which is converted to the dimensionless form 
%   u_t = gamma*f(u,v) + del^2 u 
%   v_t = gamma*g(u,v) + d*del^2 v 
% with u(r,0) and v(r,0) fixed and no-flux boundary conditions. 
% Solves on a 2-d "finger" domain consisting of a rectangle with a semi-
% ellipse attached to the end. 
% 
% Optionally, the user may specify a filename to save the final plot to. 
% This name should not have an extension: e.g., if you want to save the
% file as "myplot.eps", then enter filename='myplot'. 

function fingerprintSim (filename) 

% PARAMETERS %

WIDTH = 30; % width of finger 
LENGTH_RATIO = 1.6; % ratio of length of finger to width 
RECT_RATIO = .8; % ratio of length of rectangular part of finger to total
NSTEPSX = 300; % steps in y dimension 

% theta is the angle above horizontal that lines touch the boundary  
THETAMAX = .7 * pi / 2; 
THETAMIN = .4; 
YDELTA = -25; % y-coordinate of the delta(s)

gamma = 70; % domain scaling parameter 

TMAX = 25; % solve from time 0 up to Tmax 
NSTEPST = 20000; % number of time steps 

% initial conditions: u(r,0) = ic_u(r), where r = (x,y); likewise for v 
z0 = 0.1;
x0 = 0; 
y0 = -20; 
a = .15*(WIDTH/2); 
b = .1*(WIDTH/2); 
alpha = 0; % angle clockwise off vertical 
ic_u = @(x,y) max(z0, exp(-((x*cos(alpha)-y*sin(alpha)-(x0*cos(alpha)-y0*sin(alpha))).^2/(2*a^2) + ((x*sin(alpha)+y*cos(alpha)-(x0*sin(alpha)+y0*cos(alpha))).^2/(2*b^2))))) + .05*rand(size(x)); 
ic_v = @(x,y) ic_u(x,y); 

% Parameters in the Gierer-Meinhardt system 
a = .007; 
b = 1; 

% Note: Edit file diffCoeff.m to define values of d. 

% Graphical parameters 
N_UPDATE = 1; % number of times to update the plot visually 
THRESHHOLD = [40, 5]; % see plot_finger.m 

%%%%%%%%%%


% Check parameter values 
if RECT_RATIO >= 1 || RECT_RATIO <= 0
    error('RECT_RATIO must be between 0 and 1.') 
end
if any(THRESHHOLD < 0) || THRESHHOLD(1) + THRESHHOLD(2) > 99
	error('Invalid threshhold values.');
end

dt = TMAX / NSTEPST; % time step 
h = WIDTH / NSTEPSX; % grid step size 

% Set steps in each direction, possibly slightly altering length and 
% RECT_LENGTH, so as to keep the step size h constant. 
NSTEPSY = round(LENGTH_RATIO * WIDTH / h); 
LENGTH = h * NSTEPSY; 
NSTEPSRECT = round(NSTEPSY * RECT_RATIO); 
RECT_LENGTH = h * NSTEPSRECT; 

% Create domain. C = (u,v); so C(:,:,1) = u, C(:,:,2) = v. Further, 
% Make sure theta is in appropriate range 
THETAMAX = min(THETAMAX, .9 * pi / 2); 
THETAMAX = max(THETAMAX, 0); 
THETAMIN = max(THETAMIN, 0);
THETAMIN = min(THETAMIN, THETAMAX); 
if THETAMAX > .75 * pi / 2 
	warning('MATLAB:fingerprintSim', ['Values of theta near pi/2 '...
		'may cause numerical instability.']);
end

% C(x,y,U) = the value of u at point (x,y). 
% The domain is indexed U(x,y) where y runs from the base of the finger
% segment to the tip, and x runs from left to right as seen when facing the
% fingerprint. That is, U(1,1) is the point at the bottom left corner of
% the finger segment, and U(size(U,2)/2,size(U,1)) is the tip of the 
% finger. 
U = 1; V = 2; 
[X, Y] = meshgrid (linspace(-WIDTH/2, WIDTH/2, NSTEPSX+1), ...
                linspace(-RECT_LENGTH, LENGTH-RECT_LENGTH, NSTEPSY+1)); 
X = X'; Y = Y'; 
[YMAX, XMIN, XMAX, OUTSIDE] = fingerBoundary(X, Y, RECT_LENGTH, THETAMAX); 

% Convert ydelta to integer index  
YDELTA = min(YDELTA, Y(1,end-1)); 
IYDELTA = find(Y(1,:) >= YDELTA, 1);    
IYDELTA = max([IYDELTA, 2]); 
IYDELTA = min([IYDELTA, YMAX(1) - 1, YMAX(end) - 1]); 
% Interpolate theta linearly, copy for both sides  
ysplit = min(YMAX(1), YMAX(end)); 
theta = [zeros(1,IYDELTA-2), linspace(THETAMIN, THETAMAX, ysplit-IYDELTA)]; 
theta = [theta; theta]; 

% Diffusion coefficient d 
[d_funcs, d_nsteps] = diffCoeff(TMAX, NSTEPST, WIDTH/2, min(min(Y)), ...
						max(max(Y)), IYDELTA); 
% Check CFL condition for stability 
dcheck (d_funcs, X, Y, dt/h^2);

% Implement initial conditions 
C(:,:,U) = ic_u(X,Y); 
C(:,:,V) = ic_v(X,Y); 
C = max(C, 10^-3); % avoid zeros on first pass through "update" 
C(OUTSIDE) = 0; 

% Plot initial 
figure 
plot_finger(C(:,:,U), X, Y, YMAX, OUTSIDE, 0, THRESHHOLD);
hold off 
time = 0; 

for id = 1:numel(d_funcs) 
	
	d = d_funcs{id}(X, Y); 
	nstepst = d_nsteps(id); 
	
	nupdate = ceil(N_UPDATE * nstepst / NSTEPST);
	nupdate = min(nupdate, nstepst); 
	nupdate = max(nupdate, 1); 
	nt = round(nstepst / nupdate); 
	nupdate = floor(nstepst / nt); 
	
	for it = 1:nupdate
	
		if it == nupdate 
			nt = nstepst - nt * (nupdate - 1); 
		end
		
		C = solve (C, time, nt, dt, h, YMAX, XMIN, XMAX, a, b, d, ...
					gamma, theta); 
		
		time = time + nt * dt;
		plot_finger (C(:,:,U), X, Y, YMAX, OUTSIDE, time, THRESHHOLD);
		drawnow;  
		
	end
end

% Save image file
if nargin >= 1
	print ('-depsc', filename); 
end

end % fingerprintSim 





