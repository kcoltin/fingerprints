% dcheck.m 
% Kevin Coltin 
% 
% Checks whether the functions defining d are such that d will never be
% greater than its maximum value allowed by the stability conditions, as
% defined by the Courant number. 
% 
% Arguments: 
% d_funcs: a cell vector of functions defining d. 
% X, Y: the x- and y-coordinates of the finger domains. 
% dth2: should equal dt/h^2. 

function c = dcheck (d_funcs, X, Y, dth2) 

% Cutoff value for Courant number for ADI method with Gierer-Meinhardt 
CUTOFF = 7.5; 

dmax = 1; % default 1, since the CFL condition for the u equations
			% effectively has d = 1. 
for i = 1:numel(d_funcs) 
	dmax = max(dmax, max(max(d_funcs{i}(X,Y)))); 
end

c = 4 * dmax * dth2; 

if c >= CUTOFF
	error(['Scheme will be numerically unstable: Courant number ', ...
		'4*d*dt/h^2 should be less than %.2g.\nCourant number = ', ...
		'%.3g'], CUTOFF, c); 
end

end










