% solve.m 
% Kevin Coltin 
%
% Solves a few time steps of the reaction-diffusion PDES governing 
% fingerprint formation. This function serves an an interface between
% fingerprintSim and mxFDupdate/FDupdate.mex*: its main function is to 
% check the correct size of arguments, to avoid memory access errors. 

function Cout = solve (C, time, nt, dt, h, ymax, xmin, xmax, a, b, d, gamma, theta) 		

% Check arguments. Note that this only checks the size of arguments, 
% not their values - checking values is done in Fortran. Checking size
% must be done here to avoid memory errors in the Mex gateway function. 
		
% C must be n x m x 2 
dims = size(C); 
if numel(dims) ~= 3 
	error('C must be three dimensional.'); 
elseif dims(1) < 5 || dims(2) < 5
	error('First two dimensions of C are too small.'); 
elseif dims(3) ~= 2 
	error('Third dimension of C must equal 2.'); 
end 

% Check dimensions of domain boundary parameters 
% Note: The values of the three vectors (e.g. if ymax is compatible 
% with xmin and xmax) is not checked here, since it is guaranteed by
% fingerBoundary and would be complicated to check again.
ymax = ymax(:); xmin = xmin(:); xmax = xmax(:); % coerce to column vectors 
if any([ymax; xmin; xmax] ~= mod([ymax; xmin; xmax], 0)) ...
	|| any([ymax; xmin; xmax] <= 0) 
	error('ymax, xmin, and xmax must be integers.'); 
elseif numel(ymax) ~= dims(1) 
	error('Dimension of ymax must equal first dimension of C.'); 
elseif numel(xmin) ~= dims(2) || numel(xmax) ~= dims(2) 
	error('Dimension of xmin and xmax must equal second dimension of C.'); 
end 

% Check d 
if ndims(d) ~= 2 
	error('d must be two dimensional.');
elseif any(size(d) ~= dims(1:2)) 
	error('d must be the same size as the first two dimensions of C.'); 
end 

% Check theta 
if ndims(theta) ~= 2 
	error('theta must be two dimensional.'); 
elseif size(theta, 1) ~= 2 
	error('First dimension of theta must equal 2.'); 
elseif size(theta, 2) ~= min(ymax(1), ymax(end)) - 2
	error('Second dimension of theta must equal min(ymax(1),ymax(end))-1.');
end 

Cout = FDupdate(C, time, nt, dt, h, ymax, xmin, xmax, a, b, d,	...
				gamma, theta);
					
end 
					
