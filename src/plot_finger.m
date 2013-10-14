% plot_finger.m 
% Kevin Coltin 
% 
% Plots a "fingerprint" in formation, by plotting the concentrations of the
% primary morphogen A. 
% 
% The concentrations are plotted on a relative scale, obtained by dividing
% A / max(max(A)). This is clearly not very scientific, so we should 
% consider changing it eventually. 
% 
% Arguments: 
% A: Matrix of morphogen concentrations. 
% X, Y: Matrices of size(A) for the x- and y-coordinates on the domain. 
% YMAX: Boundary of finger domain 
% OUTSIDE: Indices of points outside the domain. 
% time: Time in the simulation, for including on the plot legend. 
% threshhold: If threshhold = [x,y], then the largest y% of values and
%	the smallest x% will be ignored when setting the bounds by normalizing
%	the values on [0,1]. That is, the values will be normalized so that
%	[a,b] maps to [0,1], where a is the xth percentile value in A and b is 
%	the (100-y)th percentile value. 

function plot_finger (A, X, Y, YMAX, OUTSIDE, time, threshhold) 

if nargin < 7
	threshhold = [2.5, 2.5]; 
elseif numel(threshhold) ~= 2
	error('Threshhold must be a vector of length 2.'); 
end

% Convert to zeros outside domain 
A(OUTSIDE) = 0; 

% Threshhold cutoff 
if threshhold(1) > 0
	minval = prctile(reshape(A, numel(A), 1), threshhold(1)); 
	A(A < minval) = minval; 
end
if threshhold(2) > 0
	maxval = prctile(reshape(A, numel(A), 1), 100 - threshhold(2)); 
	A(A > maxval) = maxval; 
end

% Create image: normalize concentrations on [0, 1]. Reverse them so that 
% 0's show up as white and 1's as black. 
if any(any(A)) ~= 0
	img = A - min(min(A)); 
	img = 1 - img / max(max(img)); 
else
	img = ones(size(A)); 
end

% Reensure it's all zeros outside domain 
img(OUTSIDE) = 1; 

img3 = ones([size(img) + 2, 3]); % add 2 to make room for border 
for i = 1:3
	img3(2:end-1,2:end-1,i) = img; 
end


% Add lines delineating domain boundary 
midpt = find(YMAX == max(YMAX), 1); 
for i = 1:midpt
	img3(i+1,YMAX(i)+2:YMAX(i+1)+2,1) = 0; 
end
for i = midpt+1:size(A,1)
	img3(i+1,YMAX(i)+2:YMAX(i-1)+2,1) = 0; 
end
img3(1,1:YMAX(1)+2,1) = 0; 
img3(:,1,1) = 0; 
img3(end,1:YMAX(end)+2,1) = 0; 

img3 = imrotate(img3, 90); 
image([2*X(1,1) - X(2,1); X(:,1); 2 * X(end,1) - X(end-1,1)],...
	[2*Y(1,end) - Y(1,end-1), fliplr(Y(1,:)), 2*Y(1,1) - Y(1,2)], img3); 
axis equal 

set(gca,'YDir', 'normal');
title({'Fingerprint formation: Morphogen concentrations', ...
		sprintf('Time = %.2f', time)}); 

	
end

