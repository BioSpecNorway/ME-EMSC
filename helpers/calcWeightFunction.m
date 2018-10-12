function [weights] = calcWeightFunction(wn, options)
    % Build weights for the reference spectrum, implemented as hyperbolic tangens functions 
    % alternating between values 0 and 1. Choose between 3 or 4 inflection points. 
    % 
    %  ---------------------------------------------------------------------------------
    %  Written by: 
    %  Johanne Solheim, Tatiana Konevskikh, Achim Kohler                                             
    %                                                                                               
    %  Faculty of Science and Technology (REALTEK)                                   
    %  Norwegian Unversity of Life Sciences (www.nmbu.no)                                                                                                                             
    %                                                                                               
    %  Post address:                                                                                                                                                                         
    %  PO Box 5003, 1432 Aas, Norway                                                               
    %                                                                                               
    %  First version: 2016                                                                      
    %  Revision: 2018  
    %  ---------------------------------------------------------------------------------
    %
    %  Input:
    %  wn       - Wavenumbers corresponding to spectra in correction (coloumn vector)  
    %  options  - Options for correcton (struct), containing fields Weights_InflectionPoints and Weights_Kappa:
    %             Weights_InflectionPoints: turning points in descending order (array type double), each pair corresponds to a chemically
    %             active area (weighted to one) 
    %             Weights_Kappa: listing steepness of function at corresponding inflection point (array type double) 
    %             Setting Weights_Kappa to 0 gives straight edges, with no smoothing. 
    %  
    %  Output: 
    %  weights  - Weights for applying on reference spectrum in ME-EMSC correction (row vector)
    
    infl1=options.Weights_InflectionPoints{1,1}(1); 
    infl2=options.Weights_InflectionPoints{1,1}(2);
    infl3=options.Weights_InflectionPoints{1,2}(1);
    infl4=options.Weights_InflectionPoints{1,2}(2); 
    
    if infl1<infl3
        msg = 'Error in weight function: Weights_InflectionPointsbers should be given in decreasing order!';
        error(msg)
    end 
    
    
    Kappa1 = options.Weights_Kappa{1,1}(1);
    Kappa2 = options.Weights_Kappa{1,1}(2);
    Kappa3 = options.Weights_Kappa{1,2}(1);
    Kappa4 = options.Weights_Kappa{1,2}(2); 

    lim = 0.0015;
 
    % Set up hyperbolic tanget function 
    extensionN = 15;
    delta_step = 0.094; 
    x=-extensionN:delta_step:extensionN;
    num_points=size(x,2);
    
    hyptan=@(a,x_range) 0.5*(tanh(a*x_range)+1); 

    % Extend wavenumber region 
    dwn = wn(2)-wn(1); 
    wn_tmp = [((wn(1)-num_points*dwn):dwn:(wn(1)-dwn))'; wn; ((wn(end)+dwn):dwn:(wn(end)+num_points*dwn))']; 
    
 
     % Find index of inflection points 
    [~,i1]=min(abs(wn_tmp-infl1));
    [~,i2]=min(abs(wn_tmp-infl2));
    [~,i3]=min(abs(wn_tmp-infl3));
    

    x1 = x; 
    x2 = x; 
    x3 = x; 
    
    if not(infl4) 
        % First patch: ones 
        p1 = ones(1, length(wn_tmp(1:i3-floor(num_points/2)))); 
        
        % Find index of end of patch 2, check if overlapping with start of patch 4
        i_end_p2 = length(p1) + num_points; 
        i_start_p4 = i2 - floor(num_points/2); 
        if i_end_p2 >= i_start_p4
            len_between = i2-i3; % Find number of points between the two inflection points 
            new_extension = x(floor(num_points/2)+floor(len_between/2)); % Find the new extension of the patch 
            x1 = -extensionN:delta_step:new_extension; % Make new x for patch 2 
            x2 = -new_extension:delta_step:extensionN; % Make new x for patch 4 
        end 
        
        % Second patch: hyperbolic tangent with infl3 and Kappa3
        p2 = -hyptan(Kappa3, x1)+1; 

        
        if p2(end)>lim % LEGG TIL 
            msg = ('Weight function: inflection point no. 2 and 3 are too close.'); 
            warning(msg)
        end 
        
        % Third patch: zeros 
        if i_start_p4-i_end_p2>=0
            p3 = zeros(1,i_start_p4-i_end_p2);
        else 
            p3 = [];
        end 
        
        i_end_p4 = i2 + floor(num_points/2); 
        i_start_p6 = i1 - floor(num_points/2); 
        if i_end_p4 >= i_start_p6
            len_between = i1-i2; % Find number of points between the two inflection points 
            %len_between = 237; % REMOVE THIS!!!
            
            new_extension = x(floor(num_points/2)+floor(len_between/2)); % Find the new extension of the patch 
            x2 = x2(1):delta_step:new_extension; % Make new x for patch 1 
            x3 = -new_extension:delta_step:extensionN; % Make new x for patch 2 
        end 
        
        % Fourth patch: hyperbolic tangent with infl2 and Kappa2
        p4 = hyptan(Kappa2, x2); 
        
        if p4(1)>lim 
            msg = ('Weight function: inflection point no. 2 and 3 are too close.');
            warning(msg)
        elseif 1-p4(end)>lim 
            msg = ('Weight function: inflection point no. 1 and 2 are too close.');
            warning(msg)
        end 
                
        % Fifth patch: ones 
        if i_start_p6-i_end_p4>=0
            p5 = ones(1,i_start_p6-i_end_p4);
        else 
            p5 = [];
        end 
        
        % Sixth patch: hyperbolic tangent with infl1 and Kappa1
        p6 = -hyptan(Kappa1, x3)+1;
        
        % Seventh patch: zeros 
        p7 = zeros(1, length(wn_tmp)-ceil(num_points/2)-i1); 
        
        % Weights 
        weights = [p1 p2 p3 p4 p5 p6 p7]; 
        
        diff = length(weights)-length(wn_tmp); 
        if diff>0
            weights = weights(1:end-diff);
        elseif diff<0
            weights = [weights repmat(weights(end), [1 -diff])]; 
        end 
    else
        [~,i4]=min(abs(wn_tmp-infl4)); 
        x4 = x; 
        % First patch: zeros
        p1 = zeros(1, length(wn_tmp(1:i4-floor(num_points/2)))); 
        
        % Find index of end of patch 2, check if overlapping with start of patch 4
        i_end_p2 = length(p1) + num_points; 
        i_start_p4 = i3 - floor(num_points/2); 
        if i_end_p2 >= i_start_p4
            len_between = i3-i4; % Find number of points between the two inflection points 
            %len_between = 237; % REMOVE THIS!!!
            new_extension = x(floor(num_points/2)+floor(len_between/2)); % Find the new extension of the patch 
            x1 = -extensionN:delta_step:new_extension; % Make new x for patch 1 
            x2 = -new_extension:delta_step:extensionN; % Make new x for patch 2 
        end 
        
        % Second patch: hyperbolic tangent with infl4 and Kappa4
         p2 = hyptan(Kappa4, x1); 
        
        if 1-p2(end)>lim 
            msg = ('Weight function: inflection point no. 3 and 4 are too close.');
            warning(msg)
        end 
         
        % Third patch: ones 
        if i_start_p4-i_end_p2>=0
            p3 = ones(1,i_start_p4-i_end_p2);
        else 
            p3 = [];
        end 
        
        % Find index of end of patch 4, check if overlapping with start of patch 6
        i_end_p4 = i3 + floor(num_points/2); 
        i_start_p6 = i2 - floor(num_points/2); 
        if i_end_p4 >= i_start_p6
            len_between = i2-i3; % Find number of points between the two inflection points 
            new_extension = x(floor(num_points/2)+floor(len_between/2)); % Find the new extension of the patch 
            x2 = x2(1):delta_step:new_extension; % Make new x for patch 1 
            x3 = -new_extension:delta_step:extensionN; % Make new x for patch 2 
        end 
        
        % Fourth patch: hyperbolic tangent with infl2 and Kappa2
        p4 = -hyptan(Kappa3, x2)+1; 
        
        if 1-p4(1)>lim 
            msg = ('Weight function: inflection point no. 3 and 4 are too close.');
            warning(msg)
        elseif p4(end)>lim 
            msg = ('Weight function: inflection point no. 2 and 3 are too close.');
            warning(msg)
        end 
        
        % Fift patch: zeros 
        if i_start_p6-i_end_p4>=0
            p5 = zeros(1,i_start_p6-i_end_p4);
        else 
            p5 = [];
        end 
        
        % Find index of end of patch 6, check if overlapping with start of patch 8
        i_end_p6 = i2 + floor(num_points/2); 
        i_start_p8 = i1 - floor(num_points/2); 
        if i_end_p6 >= i_start_p8
            len_between = i1-i2; % Find number of points between the two inflection points 
            new_extension = x(floor(num_points/2)+floor(len_between/2)); % Find the new extension of the patch 
            x3 = x3(1):delta_step:new_extension; % Make new x for patch 1 
            x4 = -new_extension:delta_step:extensionN; % Make new x for patch 2 
        end 
        
        % Sixth patch: hyperbolic tangent with infl3 and Kappa3
        p6 = hyptan(Kappa2, x3); 
        
        if p6(1)>lim 
            msg = ('Weight function: inflection point no. 2 and 3 are too close.'); 
            warning(msg)
        elseif 1-p6(end)>lim 
            msg = ('Weight function: inflection point no. 1 and 2 are too close.'); 
            warning(msg)
        end 
        
        % Seventh patch: ones
        if i_start_p8-i_end_p6>=0
            p7 = ones(1,i_start_p8-i_end_p6);
        else 
            p7 = [];
        end 
        
        % Eight patch: hyperbolic tangent with infl1 and Kappa1
        p8 = -hyptan(Kappa1, x4)+1;
        
        % Nineth patch: zeros
        p9 = zeros(1, length(wn_tmp)-ceil(num_points/2)-i1); 
        
        % Weights 
        weights = [p1 p2 p3 p4 p5 p6 p7 p8 p9]; 
        
        diff = length(weights)-length(wn_tmp); 
        if diff>0
            weights = weights(1:end-diff);
        elseif diff<0
            weights = [weights repmat(weights(end), [1 -diff])]; 
        end 
    end 
    
    weights = weights((num_points+1):end-(num_points));
end 
    

    
 