function [in, current_element] = localFind(current_point, current_element)
    
    in =  triIntersect(current_element.Faces, current_point);
    if in
        return
    end
    
    [eff_list] = efficient_list(current_element,current_point,current_element.sphere_radius);
                        
    num_elements = length(eff_list);
    j=1;
    while ~in && j<=num_elements
        current_element = eff_list(j);
        in =  triIntersect(current_element.Faces, current_point);
        j = j+1;               
    end
end