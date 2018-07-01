function [in, Element] = globalFind(PartArr, seedElement, current_point)
%% Globally Locate a Point in an Element
% This function locates a point globally. It utilises the fact that the
% path step length should be less than the maximum distance between any two
% nodes in any element. This should be obvious as, to get any detail in the
% paths, it is necessary to see hw the path moves through elements and is
% influenced by elements surrounding it.
% The algorithm excludes any nodes that arent in the sphere of influence.
% The elements that the remaining nodes belong to are then tested.
% A future update of this will eleminate entire parts at a time by taking
% the convext hull of a set of elements and ray casting the boundary.

%The stuck test was used to get around some edge case that kept causing the
%program to hang. I cant remember what tat was at this stage.

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
        eff_list = efficient_list(seedElement,current_point,initial_radius);
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
                in = triIntersect(Element.Faces, current_point');
                k = k+1;
            end
        end
    end
