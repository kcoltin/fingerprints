% fingerBoundary.m 
% Kevin Coltin 
% 
% Returns the indices defining the boundary of a finger domain. 
% ymax, the returned value, is a 1 x (NSTEPSX+1) vector whose ith entry is
% equal to max{j : (i,j) is inside the finger domain}.  
% 
% Also returns xmin and xmax, the minimium and maximum x coordinates within
% the domain for each given y. 
% 
% Note: 
% -This function guarantees that ymax has no unique maximum. That is, the
% two (or more) central points of the fingertip are of equal height. This
% prevents complications with the iyL and iyR indices in 
% central_difference.m. It is also done because otherwise the no-flux
% boundary conditions would not be well defined at such points. 
% -The function guarantees that max(ymax) == size(X,2) == size(Y,2). 

function [ymax,xmin,xmax,outside] = fingerBoundary(X, Y, rectLength, theta)

x = X(:,1); 
y = Y(1,:); 

% Force theta to be within acceptable bounds. 
% This also checks that theta is entered as a scalar- it should be the max
% value of theta, in order to determine the angle at which the semiellipse
% touches the rectangle, not the whole 2 x n matrix of theta.  
if numel(theta) ~= 1 
	error('theta must be a scalar.'); 
end
THETA_MAX = .9 * pi / 2; % almost pi / 2
theta = min(theta, THETA_MAX);
theta = max(theta, 0); 

% Equation of ellipse is x^2/a^2 + y^2/b^2 = 1. 
w = 2 * x(end); % width of domain 
h = y(end); % height of tip segment of finger 

% solve a quadratic equation for b 
b = h*(8*h^2-w^2*tan(theta)^2+2*h*w*tan(theta))/(16*h^2-w^2*tan(theta)^2);
a = b*cot(theta)*sqrt((2*b*h-h^2)/(h^2-2*b*h+b^2));
y0 = h - b;

% Construct ellipse segment based on values a, b, y0: 
% y-coordinates of the boundary 
ymaxval = y0 + b / a * sqrt((a + x) .* (a - x)); %y= y0 + b*sqrt(1-x^2/a^2)

% Begin adjusting boundary to ensure certain conditions are met: 

% This is to correct for if ymaxval(i) is slightly *greater* than the
% maximum value of y or *less* than zero due to roundoff. 
ymaxval = min(ymaxval, max(y)); 
ymaxval = max(ymaxval, 0); 

% Find the index of the value in Y closest to the coordinates ymaxval 
ymax = zeros(size(x)); 
for i = 1:numel(x)
	
	[~, index] = histc(ymaxval(i), y); 
	
	if index < numel(y) && y(index+1) - ymaxval(i) < ymaxval(i) - y(index) 
		index = index + 1; 
	end
	ymax(i) = index; 
	
end

% ensure the tip is at least three wide (this must be done so that boundary
% conditions can be implemented implicitly in fortran code) 
imax = find(ymax == max(ymax)); % indices of maximum values of y 
if numel(imax) == 1 
	ymax([imax-1, imax+1]) = max(ymax); 
elseif numel(imax) == 2 
	% Add another point to the top row, on the side that will make it more
	% symmetric 
	if imax(2) - numel(ymax) / 2 > numel(ymax) / 2 - imax(1)
		ymax(imax(1)-1) = max(ymax); 
	else
		ymax(imax(2)+1) = max(ymax); 
	end
end

% Ensure that max(ymax) is equal to size(Y,2): that is, the fingertip
% reaches all the way to the top of the domain. This is important for the
% Fortran functions. 
diff = size(Y,2) - max(ymax); 
ymax = ymax + diff; 

% Ensure that ymax is at least 2 at each point. This is necessary for
% Fortran function val_bottom_nbc; it is also necessary for other functions
% that ymax be at least 1 at all points. 
ymax = max(ymax, 2); 

% Create indices xmin and xmax. Note: the transpose is necessary to make
% them column vectors, which is necessary for LUdecomp to work correctly. 
xmin = zeros(size(y))'; 
xmax = zeros(size(y))'; 

for iy = 1:numel(y) 
	xmin(iy) = find(ymax >= iy, 1); 
	xmax(iy) = find(ymax >= iy, 1, 'last'); 
end

if any(ymax < 1) 
	error('rectLength is too small.')
end

if ymax(1) == max(ymax) % Note 1 
	NSTEPSX = numel(x) - 1; 
	NSTEPSY = numel(y) - 1; 
	ymax([1,end]) = NSTEPSY; 
	xmin(end) = 2; 
	xmax(end) = NSTEPSX; 
end

outside = getOutsideIndices (ymax, xmin, xmax);

end 


% Notes 
% 1. We need to make the top row not all equal: this prevents errors that 
% arise in the Fortran code if YMAX is all	exactly the same. 






