function [efficient_list] = efficient_list(test_list,...
                            current_point,radius_of_influence)
%% Efficient List
% This function uses the sphere of influence to eliminate elements outside
% a useful boundary. The radius of is boundary is calculated by taking the
% maximal distance between two nodes in the current element and halving
% it.
    
    x = current_point(1);
    y = current_point(2);
    z = current_point(3);
    local_nodes = [test_list(:).nodes];

    [point_to_node_dist,sort_idx] = sort([sqrt((x-[local_nodes.xCoordinate]).^2 + ...
    (y-[local_nodes.yCoordinate]).^2 + (z-[local_nodes.zCoordinate]).^2)]);

    local_nodes = local_nodes(sort_idx(point_to_node_dist<=radius_of_influence));
    
    num_nodes = length(local_nodes);
    if ~num_nodes
        efficient_list = [];
        return
    end
    
    surrounding_elements(num_nodes).elements = [];
    
    for k = 1:num_nodes
        efficient_nodal_els = [local_nodes(k).inelement(:)]';
        surrounding_elements(k).elements = efficient_nodal_els([efficient_nodal_els(:).importeddata]>0);
    end
    
    efficient_list = setdiff([surrounding_elements(:).elements], test_list,'stable');
