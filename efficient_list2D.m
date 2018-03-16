function [efficient_list] = efficient_list2D(test_list,...
                            current_point,radius_of_influence)
    
    x = current_point(1);
    y = current_point(2);

    local_nodes = [test_list(:).nodes];

    [point_to_node_dist,sort_idx] = sort([sqrt((x-[local_nodes.xCoordinate]).^2 + ...
    (y-[local_nodes.yCoordinate]).^2)]);

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