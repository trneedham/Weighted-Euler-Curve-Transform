function SWECT = complex_to_weighted_ECT(complex,num_directions,num_steps,method,window,normalization_method)

% Inputs: 
% complex = weighted simplicial complex representing a 2D greyscale image.
%           This is a structure containing
%           complex.V = vertex positions
%           complex.E = edge connection data
%           complex.F = triangle data
%           complex.V_weights, complex.E_weights, complex.F_weights =
%           weights on the above
%
%           Such a complex can be built from an image via
%           'build_weighted_complex.m'
%
% num_directions = number of directions to compute Euler curves for
% num_steps = number of steps in Euler curve
% method = smoothing method, e.g. 'gaussian', 'movmean', etc. Any method
%          valid for the 'smoothdata' Matlab function will work here
% window = smoothing window for 'smoothdata'
% normalization_method = method to normalize the weights on the complex.
% The options here are:
%           'none' = leave the weights as-is
%           'max' = normalize every value by the max of all values
%           'total' = normalize everything by the sum of all weights; i.e., vertex
%                     weights are all divided by the sum of all vertex weights
%           'ECT' = set all weights equal to one. This is just the (S)ECT

% Output: num_steps x num_directions matrix SWECT. Each SWECT(:,j) is a smoothed 
% weighted Euler curve for the jth direction.

if nargin == 5
    normalization_method = 'none';
end

% Create directions

rotstep=num_directions;
theta=-pi:2*pi/rotstep:pi;
d1=cos(theta);d2=sin(theta);
d=[d1;d2];

% Extract simplicial complex data

V = complex.V;
E = complex.E;
F = complex.F;

if strcmp(normalization_method,'none')
    V_weights = complex.V_weights;
    E_weights = complex.E_weights;
    F_weights = complex.F_weights;
elseif strcmp(normalization_method,'max')
    V_weights = complex.V_weights/max(complex.V_weights);
    E_weights = complex.E_weights/max(complex.V_weights);
    F_weights = complex.F_weights/max(complex.V_weights);
elseif strcmp(normalization_method,'total')
    V_weights = complex.V_weights/sum(complex.V_weights);
    E_weights = complex.E_weights/sum(complex.E_weights);
    F_weights = complex.F_weights/sum(complex.F_weights);
elseif strcmp(normalization_method,'ECT')
    V_weights = ones(length(complex.V),1);
    E_weights = ones(length(complex.E),1);
    F_weights = ones(length(complex.F),1);
end

% Center the simplicial complex at the origin. 

Z = V;
Z=Z-repmat(mean(Z),length(Z),1);

% Normalize the complex so that its max vertex distance from the origin is
% 1. Remove these lines to use unnormalized complexes.
% If these lines are removed, then the 'weighted_euler_curve.m' function
% will also have to be changed, since it currently assumes normalized
% complexes!

r = max(vecnorm(Z'));
Z = Z/r;
    

% Create the Euler curve transform by computing an Euler curve in each
% direction

SWECT = zeros(num_steps,num_directions);

for i=1:length(d)-1
    fun=Z*d(:,i);
    EC=weighted_euler_curve(Z,E,F,V_weights,E_weights,F_weights,fun,num_steps);
    if strcmp(method,'none')
        SWECT(:,i) = EC(1:num_steps,2);
    else
        SWECT(:,i) = smoothdata(EC(1:num_steps,2),method,window);
    end
end