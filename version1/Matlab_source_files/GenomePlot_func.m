% GenomePlot_func
% Ver 2.00 - switch figures 1 and 2 regarding isochoric domains.
% Ver 3.00 - version 2 was used to do the plotting for the paper and it works fine. Here i am trying to create a multi-chromosomal display.
% Ver 4.00 - plotting all lines together as a matrix to reduce complexity of the plot, which interfered with Linux
% Website: http://code.google.com/p/isoplotter/

function GenomePlot_func(plot_type, iso_size, min_size, M2, Y_data, curr_chr_len_all_prop)

    %Get the size, GC%, homogeneity indication, and isochoric indication for each domain
    x = M2(M2(:,4)>=min_size,4)'; %size
    y = M2(M2(:,4)>=min_size,5)'; %GC%
    hom_domain = M2(M2(:,4)>min_size,7)'; %homogeneous
    iso_domain = x>=iso_size; 

%     disp(y(find(iso_domain))) %display isochores GC%
%     disp([sum(y(iso_domain)<0.2) sum(y(iso_domain)>0.6)]);    %are there isochoric domains with extreme values

    color_range = 100;

    %Colorful plot
    if plot_type==1
        color_vec = [1 0 0]; %Set isochoric domains color 
    
    %Homogeneous/non homogeneous
    elseif plot_type==2

        fractions =   [0.2 0.2 0.2 0.2 0.2];
        color_vec_a = [ones(fractions(1)*color_range,1)      ones(fractions(1)*color_range,1)       zeros(fractions(1)*color_range,1)];     %yellow [1 1 0]
        color_vec_b = [zeros(fractions(2)*color_range,1)     0.5.*ones(fractions(2)*color_range,1)  ones(fractions(2)*color_range,1)];      %light blue [0 0.5 1]
        color_vec_c = [ones(fractions(3)*color_range,1)      zeros(fractions(3)*color_range,1)      zeros(fractions(3)*color_range,1)];     %red [0 1 0]
        color_vec_d = [zeros(fractions(4)*color_range,1)     ones(fractions(4)*color_range,1)       zeros(fractions(4)*color_range,1)];     %green [0 1 0]
%         color_vec_e = [0.5.*ones(fractions(5)*color_range,1) 0.5.*ones(fractions(5)*color_range,1)  0.5.*ones(fractions(5)*color_range,1)]; %gray [0 0 1]
        color_vec_e = [ones(fractions(5)*color_range,1) 0.5.*ones(fractions(5)*color_range,1)  zeros(fractions(5)*color_range,1)]; %orange [1 0.5 0]
        
        color_vec = [color_vec_a' color_vec_b' color_vec_c' color_vec_d' color_vec_e' ]';

        hom_domain(iso_domain==0)=0; %Take all non isochoric domains and mark them as non homogeneous to set their color to white
        
    %Isochoric/non isochoric
    elseif plot_type==3
        color_vec = [1 0 0];
        hom_domain(iso_domain==0)=0; %Take all non isochoric domains and mark them as non homogeneous to set their color to white
    end;
        
%     colormap(color_vec);
%     colorbar

    %QC - make sure all fractions are there
    if exist('fractions','var')
        if sum(fractions)~=1
            disp('Fractions do not sum to 1. returning.');
            return;
        end;
    end;

    %Print chr size
    x_cum = cumsum(x);
    x_rel = cumsum([0 x./max(x_cum)]);
    x_rel = x_rel * curr_chr_len_all_prop;% new in version
%     disp(max(x_cum));

    %Print color domains
    final_color = [];
    for i=1:length(x) 
        curr_color = round(size(color_vec,1)*y(i)); %Multiply the color in GC%
        
        %Make sure we stay in range
        if curr_color==0 
            curr_color = 1;
        elseif curr_color>size(color_vec,1)
            curr_color = size(color_vec,1);
        end;
    
        %Homogeneous
        if hom_domain(i)
            %Check if > iso size
            if iso_domain(i)
                %Long Isochoric
%                 patch('XData',[x_rel(i); x_rel(i+1); x_rel(i+1); x_rel(i)], 'YData',Y_data, 'FaceColor', color_vec(curr_color,:), 'EdgeColor', color_vec(curr_color,:)); %Long iso domains - colorful (ver 3.00)
                final_color(i,:) = color_vec(curr_color,:);
    %             disp([i x_rel(i) x_rel(i+1) x(i) y(i)]);
            %Short Isochoric
            else
%                 patch('XData',[x_rel(i); x_rel(i+1); x_rel(i+1); x_rel(i)], 'YData',Y_data, 'FaceColor', 'k', 'EdgeColor', 'k'); % Short are black (ver 3.00)
                final_color(i,:) = [0 0 0];
            end;
        else
            %Non-homogeneous
%             patch('XData',[x_rel(i); x_rel(i+1); x_rel(i+1); x_rel(i)], 'YData',Y_data, 'FaceColor', 'w', 'EdgeColor', 'w' ); % Non homogeneous are white (ver 3.00)
            final_color(i,:) = [1 1 1];
        end;
    end;

    i=1:length(x); j=i+1;
    patch('XData',[x_rel(i); x_rel(j); x_rel(j); x_rel(i)], 'YData', repmat(Y_data, 1,numel(i)),...
        'FaceColor','flat',...
        'FaceVertexCData',final_color,'edgecolor','none');
end
