%Input: Finite dimensional simplicial complex representing 2D image
%       V is n1x3 vertices
%       E is n2x2 edges
%       F is n3x3 faces
%       V_weights, E_weights, F_weights are weights on V, E, F
%       fun: The function value on the vertices
%       stepsize: discrete Euler Characteristic curve into such dimensional
%       vector. (e.g. 100);
%Output: Stepsize_by_2 matrix. Euler Characteristics (second column) for scales of function values (first column).

% This code was adapted from the Matlab implementation here: https://github.com/lorinanthony/SECT
% The code there is an implementation of:
% L. Crawford, A. Monod, A.X. Chen, S. Mukherjee, and R. Rabadán. 
% Predicting clinical outcomes in glioblastoma: an application of topological and functional data analysis. 
% Journal of the American Statistical Association. In Press.

function curve = weighted_euler_curve_tumors(V,E,F,V_weights,E_weights,F_weights,fun,stepsize)

if length(V)~=length(fun)
    fprintf('The size of function should be same as the number of vertices');
    return
end

fe=zeros(size(E,1),1);
ff=zeros(size(F,1),1);

fe=max(fun(E)')';
ff=max(fun(F)')';

% The next step assumes that the simplicial complex has been normalized so
% that its farthest vertex from the origin is at distance one. This needs
% to be adjusted for unnormalized complexes.

threshold = -1:2/stepsize:1;

curve = zeros(stepsize+1,2);

curve(:,1)=threshold';

for i=1:length(threshold)
    vSum = sum(V_weights(find(fun<=threshold(i))));
    eSum= sum(E_weights(find(fe<=threshold(i))));
    fSum= sum(F_weights(find(ff<=threshold(i))));

    curve(i,2)=vSum-eSum+fSum;
end
