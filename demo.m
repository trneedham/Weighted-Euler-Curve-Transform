%%%%%%%%%%%%%%%%%%%%
% Demonstration of basic functionality of the WECT code.
%%%%%%%%%%%%%%%%%%%%

% This code implements the topological shape analysis technique described
% in the paper
% 'The Weighted Euler Characteristic Transform for Image and Shape
% Analysis' by Qitong Jiang, Sebastian Kurtek and Tom Needham
% See the paper for explanations of what's going on here.
% Please contact tneedham@fsu.edu with questions or to report bugs in the
% code.

%% Load greyscale images

% Add path to MNIST digit images
addpath('./Data')

% Load in MNIST Data
for j = 1:10
    MNISTdata = import_MNIST('Mnist_test.csv',j-1);
    for k = 1:100
        dataset{j,k} = MNISTdata(:,:,k)';
    end
    ['Digit ',num2str(j-1),' loaded...']
end

%% Convert the greyscale image to a weighted simplicial complex

% Pick a single image to analyze first

image = dataset{3,5};

% The function 'build_weighted_complex' returns Vertices, Edges, Faces and
% their weights. These can be stored in a struct called 'complex' to be
% passed into the next function.

[V,E,F,V_weights,E_weights,F_weights] = build_weighted_complex(image);
complex.V = V;
complex.E = E;
complex.F = F;
complex.V_weights = V_weights;
complex.E_weights = E_weights;
complex.F_weights = F_weights;

% Take a look at the result

figure(1)
clf

subplot(1,2,1)
imagesc(image)
title('Original image')
axis equal
axis off

subplot(1,2,2)
[V,E,F,V_weights,E_weights,F_weights] = build_weighted_complex(image);
triplot(F,V(:,1),V(:,2))
hold on
scatter(V(:,1),V(:,2),[],V_weights)
title('Weighted Simplicial Complex')
axis equal
axis off


%% Convert the weighted simplicial complex to its smoothed WECT

% There are several parameters to choose in the transform:

num_directions = 25; % Number of directions to filter over
num_steps = 50; % Number of time steps in each Euler curve
method = 'gaussian'; % Smoothing method for the Euler curves. 
% Any method which is valid for the Matlab function 'smoothdata' will work
% here. Use 'none' to use raw Euler curves.
window = 0.2*num_steps; % Smoothing window - larger window gives smoother curves.
normalization_method = 'none'; % Normalize the pixel values. 
% See the documentation for the 'complex_to_weighted_ECT' function for
% options.

SWECT = complex_to_weighted_ECT(complex,num_directions,num_steps,method,window,normalization_method);

% The output is a num_steps x num_directions matrix. Each SWECT(:,j) is a
% smoothed Euler curve

figure(2)
clf
plot(SWECT)
title('Smoothed Weighted Euler Curves')



%% Toy Experiment: MNIST classification

% Let's compute a pairwise distance matrix for all MNIST digits in our
% dataset. 

%% Build Complexes

% First we build simplicial complexes for each digit. This steps takes a
% while. If one is doing further analysis on these complexes, it is worth
% it to only do this once and save the result.

all_MNIST_complexes = cell(1000);

for j = 1:10
    for k = 1:100
        image = dataset{j,k};
        [V,E,F,V_weights,E_weights,F_weights] = build_weighted_complex(image);
        complex.V = V;
        complex.E = E;
        complex.F = F;
        complex.V_weights = V_weights;
        complex.E_weights = E_weights;
        complex.F_weights = F_weights;
        all_MNIST_complexes{(j-1)*100+k} = complex;
    end
    ['Digit ',num2str(j-1),' done...']
end

%% Build Weighed ECTs for each complex

% Change parameters here, if desired
num_directions = 25;
num_steps = 50; 
method = 'gaussian'; 
window = 0.2*num_steps; 
normalization_method = 'none'; 

% Build a smoothed WECT for each MNIST digit simplicial complex

SWECTs = cell(1000);

for n = 1:1000
    complex = all_MNIST_complexes{n};
    SWECTs{n} = complex_to_weighted_ECT(complex,num_directions,num_steps,method,window,normalization_method);
end

%% Compute Pairwise Euclidean Distances

% Each smoothed WECT is a num_steps x num_directions matrix. These can be
% compared with Frobenius norm

distMat = zeros(1000,1000);

for j = 1:1000
    for k = j+1:1000
        distMat(j,k) = norm(SWECTs{j}-SWECTs{k},'fro');
    end
end

distMat = distMat + distMat';

figure(3)
clf

imagesc(distMat)
title('Euclidean Distances Between Digits')

%% Compute Distances Modulo Rotation

% Since the distance computation is so fast, we can also compute a
% 'rotation invariant' version by cyclically permuting columns of one of
% the matrices and taking the min Frobenius distance. This has the effect
% of finding a registration of the image over rotations.

% For MNIST, this degrades the results since it will tend to give smaller
% distances between 6s and 9s, for example. In other applications, this
% version of distance could be quite useful.

distMat_modRotations = zeros(1000,1000);

for j = 1:1000
    for k = j+1:1000
        distMat_modRotations(j,k) = distance_RotationInvariant(SWECTs{j},SWECTs{k});
    end
end

distMat_modRotations = distMat_modRotations + distMat_modRotations';

figure(4)
clf

imagesc(distMat_modRotations)
title('Distances Mod Rotation Between Digits')

%% kNN Classification Rates

% We can compute 'Leave One Out k-Nearest Neighbor classification rates for
% these distance matrices to get a quantitative idea of how well digits are
% distinguished.

% Choose parameters
K = 1; % K in KNN
matrix = distMat_modRotations; % distance matrix to analyze

% Run the classification test
class_rates = zeros(1,10);

[~,inds] = sort(matrix,2);

neighbor_inds = inds(:,2:1+K);

neighbor_classes = ceil(neighbor_inds/100);

for j = 1:1000
    classification(j) = mode(neighbor_classes(j,:));
end

for digit = 1:10
    class_rates(digit) = sum(classification((digit-1)*100 + 1: (digit-1)*100 + 100) == digit)/100;
end

confusion = zeros(10,10);

for j = 1:10
    for k = 1:10
        confusion(j,k) = sum(classification((j-1)*100+1:(j-1)*100 + 100)==k)/100;
    end
end


overall_classification_rate = mean(class_rates)

figure(5)
clf
imagesc(confusion)
colorbar
title('Confusion Matrix')

% Note: Since the WECTs are Euclidean (if we don't use the
% rotation-invariant distance), any machine learning classifier can be run
% on the WECT signatures.
