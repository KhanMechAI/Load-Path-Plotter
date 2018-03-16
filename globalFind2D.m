function [in, Element] = globalFind2D(PartArr, seedElement, current_point)
    in = false;
    Element = seedElement;
    prev_list = seedElement;
    initial_radius = PartArr(1).maxRadius;
    tot_els = sum([PartArr(:).span]);
    stuck_test = length(prev_list);
    all_elements = [PartArr(:).elements];
    while ~in && length(prev_list) < tot_els
        [seedElement,prev_list] = surrrounding_elements(seedElement, prev_list);
        all_elements = setdiff(all_elements, prev_list);
        eff_list = efficient_list2D(seedElement,current_point,initial_radius);
        if isempty(eff_list)
            in = false;
            if stuck_test == length(prev_list)
                rand_element = randi([1,length(all_elements)],1,1);
                seedElement = all_elements(rand_element);
            end
            stuck_test = length(prev_list);
        else
            k = 1;
            while ~in && k <=length(eff_list)
                Element = eff_list(k);
                in = inpolygon(current_point(1),current_point(2), Element.verticies(:,1), Element.verticies(:,2));
                k = k+1;
            end
        end
    end