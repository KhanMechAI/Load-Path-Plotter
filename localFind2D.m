function [in, current_element] = localFind2D(current_point, current_element)
    
    in =  inpolygon(current_point(1),current_point(2), current_element.verticies(:,1), current_element.verticies(:,2));
    if in
        return
    end
    
    [eff_list] = efficient_list2D(current_element,current_point,current_element.sphere_radius);
                        
    num_elements = length(eff_list);
    j=1;
    while ~in && j<=num_elements
        current_element = eff_list(j);
        in =  inpolygon(current_point(1),current_point(2), current_element.verticies(:,1), current_element.verticies(:,2));
        j = j+1;               
    end
end