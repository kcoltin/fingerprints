% diffCoeff.m 
% Kevin Coltin 
% 
% Defines functions that set the variable d, the diffusion coefficient. 

function [d_funcs, d_nsteps] = diffCoeff (tmax, nstepst, xmax_in, ymin_in, ymax_in, iydelta_in) 

%d_funcs = {@vararch1, @vararch2, @vararch3, @homog}; 
%d_times = [0, 2, 4, 6];

d_funcs = {@homog}; 
d_times = [0];

% Parameters 
global xmax ymin ymax iydelta; % domain size parameters
xmax = xmax_in; 
ymin = ymin_in; 
ymax = ymax_in; 
iydelta = iydelta_in; 

global a b x0 y0; % ellipse parameters 
a = .3*xmax; % semi-axis in x-direction 
b = .3*xmax; % semi-axis in y-direction 
x0 = 0; % x-coordinate of center of ellipse 
y0 = -20; % y-coordinate of center of ellipse 

global DSUB DSUPER; % sub- and super-critical values of d
DSUB = 2; 
DSUPER = 6.2; 


if numel(d_funcs) ~= numel(d_times) 
	error('Length of d_funcs must equal length of d_times.');
end
if d_times(1) ~= 0 
	warning('MATLAB:diffCoeff.m', 'Changing first entry of d_times to zero.');
	d_times(1) = 0;
end

% Add ending value for last d function 
d_times(end+1) = tmax; 

if numel(d_funcs) == 1 
	d_nsteps = nstepst; 
else
	% Create number of steps for each value of d, and adjust d_times
	% accordingly
	dt = tmax / nstepst; 
	d_nsteps = round((d_times(2:end) - d_times(1:end-1)) / dt); 
	d_nsteps(end) = nstepst - sum(d_nsteps(1:end-1));
end 

end % diffCoeff 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Functions that actually define values of d 

% Arch: critical in top and bottom first, then in middle 
function d = arch (x, y) 
	global DSUB DSUPER;
	global xmax ymin ymax iydelta;
	d = zeros(size(x)); 
	ix = .7*ymin < y & y < -.3*ymax; 
	d(ix) = DSUB; 
	d(~ix) = DSUPER; 
end 

% Arch: critical in top and bottom first, gradually incrementing to the 
% middle 
function d = vararch1 (x, y)
	global DSUB DSUPER;
	global xmax ymin ymax iydelta;
	d = zeros(size(x)); 
	ix = ymin+.25*(iydelta-ymin) < y & y < .25*(ymax-iydelta); 
	d(ix) = DSUB; 
	d(~ix) = DSUPER; 
end 

function d = vararch2 (x, y)
	global DSUB DSUPER;
	global xmax ymin ymax iydelta;
	d = zeros(size(x)); 
	ix = ymin+.5*(iydelta-ymin) < y & y < .5*(ymax-iydelta); 
	d(ix) = DSUB; 
	d(~ix) = DSUPER; 
end 

function d = vararch3 (x, y)
	global DSUB DSUPER;
	global xmax ymin ymax iydelta;
	d = zeros(size(x)); 
	ix = ymin+.75*(iydelta-ymin) < y & y < .75*(ymax-iydelta); 
	d(ix) = DSUB; 
	d(~ix) = DSUPER; 
end 

% arch, d critical in middle first, then top and bottom 
function d = var2arch (x, y)  
	global DSUB DSUPER;
	global xmax ymin ymax iydelta;
	d = zeros(size(x)); 
	ix = .7*ymin < y & y < -.3*ymax; 
	d(ix) = DSUPER; 
	d(~ix) = DSUB; 
end 

% Arch: critical in top first, then in middle, then bottom 
function d = var3arch1 (x, y) 
	global DSUB DSUPER;
	global xmax ymin ymax iydelta;
	d = zeros(size(x)); 
	ix = 0 < y; 
	d(ix) = DSUPER; 
	d(~ix) = DSUB; 
end 
function d = var3arch2 (x, y)  
	global DSUB DSUPER;
	global xmax ymin ymax iydelta;
	d = zeros(size(x)); 
	ix = .7*ymin < y;
	d(ix) = DSUPER; 
	d(~ix) = DSUB; 
end 



% whorl: critical in ellipse, then everywhere
function d = whorl (x, y) 
	global a b x0 y0; 
	global DSUB DSUPER;
	global xmax ymin ymax iydelta;
	d = zeros(size(x)); 
	ix = (x - x0).^2 / a^2 + (y - y0).^2 / b^2 <= 1; 
	d(ix) = DSUPER; 
	d(~ix) = DSUB; 
end

% whorl: critical at top and bottom, then in ellipse, then everywhere 
function d = varwhorl1 (x, y) 
	global a b x0 y0; 
	global DSUB DSUPER;
	global xmax ymin ymax iydelta;
	d = zeros(size(x)); 
	ix = y < y0 - 1.2*b | y > y0 + 1.2*b; 
	d(ix) = DSUPER; 
	d(~ix) = DSUB; 
end
function d = varwhorl2 (x, y) 
	global a b x0 y0; 
	global DSUB DSUPER;
	global xmax ymin ymax iydelta;
	d = zeros(size(x)); 
	ix = (y < y0 - 1.2*b | y > y0 + 1.2*b) | ...
		(x - x0).^2 / a^2 + (y - y0).^2 / b^2 <= 1; 
	d(ix) = DSUPER; 
	d(~ix) = DSUB; 
end

% whorl: critical at everwhere except ellipse, then everywhere.  
function d = var2whorl (x, y) 
	global a b x0 y0;  
	global DSUB DSUPER;
	global xmax ymin ymax iydelta;
	d = zeros(size(x)); 
	ix = (x - x0).^2 / a^2 + (y - y0).^2 / b^2 <= .5; %*should be one, I changed it to fit the paraboloid 
	d(ix) = DSUB; 
	d(~ix) = DSUPER; 
end



% whorl: critical in ellipse, then at top and bottom, then everywhere
function d = var3whorl1 (x, y) 
	global a b x0 y0; 
	global DSUB DSUPER;
	global xmax ymin ymax iydelta;
	d = zeros(size(x)); 
	ix = (x - x0).^2 / a^2 + (y - y0).^2 / b^2 <= 1; 
	d(ix) = DSUPER; 
	d(~ix) = DSUB; 
end

function d = var3whorl2 (x, y) 
	global a b x0 y0; 
	global DSUB DSUPER;
	global xmax ymin ymax iydelta;
	d = zeros(size(x)); 
	ix = (x - x0).^2 / a^2 + (y - y0).^2 / b^2 <= 1 | y < .7*ymin | -.3*ymax < y;
	d(ix) = DSUPER; 
	d(~ix) = DSUB; 
end

% Critical everywhere except deltas, then in deltas 
function d = var4whorl (x, y)
	global DSUB DSUPER; 
	global xmax ymin ymax iydelta; 
	d = zeros(size(x)); 
	ix = abs(x) > -11 & -26 < y & y < -15;
	d(~ix) = DSUPER; 
	d(ix) = DSUB; 
end


% Loop: critical first in loop and side part of loop, then everywhere
function d = loop (x, y) 
	global a b x0 y0; 
	global DSUB DSUPER; 
	global xmax ymin ymax iydelta; 
	d = zeros(size(x)); 
	ix = (x - x0).^2 / a^2 + (y - y0).^2 / b^2 <= 1 ...
		| (x > 0 & iydelta-.5*b < y & y < iydelta+.5*b);
	d(ix) = DSUPER; 
	d(~ix) = DSUB; 
end

% Loop: critical first in loop, then top and bottom, then everywhere
function d = varloop1 (x, y) 
	global a b x0 y0; 
	global DSUB DSUPER; 
	global xmax ymin ymax iydelta; 
	d = zeros(size(x)); 
	ix = (x - x0).^2 / a^2 + (y - y0).^2 / b^2 <= 1; 
	d(ix) = DSUPER; 
	d(~ix) = DSUB; 
end

function d = varloop2 (x, y) 
	global a b x0 y0; 
	global DSUB DSUPER; 
	global xmax ymin ymax iydelta; 
	d = zeros(size(x)); 
	ix = (x - x0).^2 / a^2 + (y - y0).^2 / b^2 <= 1 | y < .6*ymin | -.5*ymax < y;
	d(ix) = DSUPER; 
	d(~ix) = DSUB; 
end

% homogeneous everywhere
function d = homog (x, ~) 
	global DSUPER;
	d = DSUPER * ones(size(x)); 
end 


%#ok<*NUSED>
%#ok<*DEFNU>





