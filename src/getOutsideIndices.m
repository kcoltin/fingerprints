% getOutsideIndices.m 
% Kevin Coltin 
% 
% Returns the indices of all points outside the finger domain. 

function outside = getOutsideIndices (ymax, xmin, xmax) 

outside = []; 
for iy = min(ymax)+1:numel(xmin)
	numbefore = (iy - 1) * numel(ymax); % number of indices before this column
	indicesonleft = numbefore + (1:xmin(iy)-1)'; 
	indicesonright = numbefore + (xmax(iy)+1:numel(ymax))'; 
	outside = [outside; indicesonleft; indicesonright];  %#ok<AGROW>
end

end 


