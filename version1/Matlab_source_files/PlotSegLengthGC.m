% PlotSegLengthGC
% Create exe file (2.5Min): tic; mcc -m -I 'D:\My Documents\Programming\MatLab\Work\' -d 'D:\Temp\' PlotSegLengthGC; toc

function PlotSegLengthGC(input_file, output_file)

    disp('Start program PlotSegLengthGC...');
    % Load file
    disp('Loading input file...');
    M = load(input_file);
    disp('File was loaded...');

    indexes = find(M(:,5)>0);    %do not plot N domains
    curr_seg = M(indexes,4);
    curr_gc = M(indexes,5);
    mean_gc = sum(curr_gc.*curr_seg)/sum(curr_seg);
    disp(['The mean GC content for the genome is ' num2str(mean_gc)]);

    %Plot the length vs GC content distribution
    figure
    plot(log(curr_seg), curr_gc, '.'); hold on;
    axis([8 16.1 0.2 0.8]);
    v = axis;
    line([v(1) v(2)],[mean_gc mean_gc], 'Color', 'r' );

    box off;
    grid off;
    ylabel('GC content');
    xlabel('Compositional domain size (log10)');

    disp(['Saving file to ' output_file]);
    print('-dtiff', output_file, '-r150');
    close all;
    disp('End program PlotSegLengthGC.');
end