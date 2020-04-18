function [V,E,F,V_weights,E_weights,F_weights] = build_weighted_complex(image)

% Input: integer-valued matrix representing a greyscale image.
% Output: simplicial complex representing the image with weights on all
% vertices, edges and faces

[height,~] = size(image);

% Find nonzero locations in the image

[row, col] = find(image);

% Define the vertex locations for the centers of the pixels. 
% We'll apply a coordinate transformation so that these are in  standard Euclidean coordinates.
% We also define the weights of these vertices by their grayscale value

V_centers = [col,height - row];

num_centers = size(V_centers,1);

V_center_weights = zeros(num_centers,1);

for j = 1:num_centers
    V_center_weights(j) = image(height - V_centers(j,2),V_centers(j,1));
end


% For each center vertex, we want to define a mini-simplicial complex like
% this:
%
%   *-------*
%   | \   / |
%       *    
%   | /   \ |
%   *-------*
%

% First define the extra vertices

V_corners = zeros(4*num_centers,2);

for j = 1:num_centers
    center_vertex = V_centers(j,:);
    V_corners((j-1)*4+1,:) = [center_vertex(1)-1/2,center_vertex(2)-1/2];
    V_corners((j-1)*4+2,:) = [center_vertex(1)+1/2,center_vertex(2)-1/2];
    V_corners((j-1)*4+3,:) = [center_vertex(1)-1/2,center_vertex(2)+1/2];
    V_corners((j-1)*4+4,:) = [center_vertex(1)+1/2,center_vertex(2)+1/2];
end

% Remove duplicates

V_corners = unique(V_corners, 'rows');

V = [V_centers; V_corners];


% Next add faces. There are four faces for each center vertex.

F = [];

for vertex_ind = 1:num_centers
    vertex = V(vertex_ind,:);
    
    NE_neighbor = vertex + [1/2,1/2];
    NE_neighbor_index = find(ismember(V,NE_neighbor,'rows'));
    NW_neighbor = vertex + [-1/2,1/2];
    NW_neighbor_index = find(ismember(V,NW_neighbor,'rows'));
    SE_neighbor = vertex + [1/2,-1/2];
    SE_neighbor_index = find(ismember(V,SE_neighbor,'rows'));
    SW_neighbor = vertex + [-1/2,-1/2];
    SW_neighbor_index = find(ismember(V,SW_neighbor,'rows'));
    
    F = [F;[vertex_ind,NE_neighbor_index,NW_neighbor_index]];
    F = [F;[vertex_ind,NW_neighbor_index,SW_neighbor_index]];
    F = [F;[vertex_ind,SW_neighbor_index,SE_neighbor_index]];
    F = [F;[vertex_ind,SE_neighbor_index,NE_neighbor_index]];
end

% Finally, extract edges from triangulation.

TO = triangulation(F,V(:,1),V(:,2));
E = edges(TO);

% Now define weights on everything. 

% Weights on the faces are easy.

for j = 1:num_centers
    F_weights((j-1)*4+1) = V_center_weights(j);
    F_weights((j-1)*4+2) = V_center_weights(j);
    F_weights((j-1)*4+3) = V_center_weights(j);
    F_weights((j-1)*4+4) = V_center_weights(j);
end

F_weights = F_weights';

% Now define weights on corner vertices by assigning highest face weight
% for containing faces

V_weights = V_center_weights;

for vertex_ind = num_centers+1:size(V,1)
    [row,~] = find(ismember(F,vertex_ind));
    V_weights = [V_weights;max(F_weights(row))];
end

% Finally define edge weights. 
    
E_weights = [];

for edge_ind = 1:size(E,1)
    edge = E(edge_ind,:);
    face_inds = find(sum(ismember(F,edge),2)==2);
    E_weights = [E_weights;max(F_weights(face_inds))];
end
