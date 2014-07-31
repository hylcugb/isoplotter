% PlotSpatialGC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    Plot Spctial GC content
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This sciprt plots the compositional domains along a chromosome
% Example: PlotSpatialGC('Z:\IsoPlotter\Examples\seg_no_ns_H.txt', 1, 'Z:\IsoPlotter\Examples\seg_no_ns_H.tif')
% 
% Create exe file (2.5Min): 
% In Windows (32bit): tic; mcc -m -I 'D:\My Documents\Programming\MatLab\Work\' -d 'D:\Temp\' PlotSpatialGC; toc
% In Windows (64bit): tic; mcc -m -I 'Z:\Tools\Matlab\Work\' -d 'Z:\IsoPlotter\' PlotSpatialGC; toc
% In Linux: tic; mcc -m -I './' -d './' PlotSegLengthGC; toc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Written by Eran Elhaik
% Written in : 11/13/12
% Ver : 2.00
% Website: http://code.google.com/p/isoplotter/
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PlotSpatialGC(input_file, chromosome_num, output_file)

    disp('Start program PlotSpatialGC...');
    
    if ~isnumeric(chromosome_num)
        chromosome_num = str2num(chromosome_num);
    end;
    
    % Load file
    disp(['Loading input file: ' input_file ' for chromosome # ' num2str(chromosome_num)]);
    M = load(input_file);
    M = M(M(:,1)==chromosome_num,:);
    disp(['File was loaded with #' num2str(size(M,1)) ' records']);

    indexes = find(M(:,5)>0);    %do not plot N domains
    curr_seg = M(indexes,4);
    curr_gc = M(indexes,5);
    mean_gc = sum(curr_gc.*curr_seg)/sum(curr_seg);
    disp(['The mean GC content for the genomic region is ' num2str(mean_gc)]);

    %Plot the length vs GC content distribution
    figure
    cx = cumsum([1 curr_seg']);
    for k=1:numel(cx)-1
        line([cx(k) cx(k+1)],[curr_gc(k) curr_gc(k)], 'Color',[0.8471 0.1608 0], 'LineWidth',5);
    end;
    hold on;
    V = axis;   %Save the vector of all axises in the graph
    line([0 V(2)],[mean_gc mean_gc], 'Color','k', 'LineWidth',1);       
    
    title(['Spatial distribution of GC content']);
    xlabel('Location along the chromosome');
    ylabel('GC content');
    box off;
    grid off;

    disp(['Saving file to ' output_file]);
    print ('-dtiff', output_file, '-r150');
    close all;
    disp('End program PlotSpatialGC.');
end

