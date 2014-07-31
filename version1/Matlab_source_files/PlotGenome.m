% PlotGenome
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    Plot Genome
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This sciprt reads a segmentation results file and plots three ideograms of the genome
% Example: PlotGenome('Z:\IsoPlotter\Examples\seg_no_ns_H.txt', 'Z:\IsoPlotter\Examples\seg_no_ns_H.tif')
% Example: PlotGenome('D:\My Documents\University\Dan Lab\01. Isochores Project\Writing\IsoPlotter\GoogleCode\IsoPlotter_2_4\Windows_executables\32bit\Example\seg_no_ns_H.txt', 'd:\seg_no_ns_H.tif')
% Example: PlotGenome('./Example/seg_no_ns_H_short.txt', './Output/PlotGenome.tif')
% 
% Create exe file (2.5Min): tic; mcc -m -I '../Tools/Matlab/Work/' -d '../IsoPlotter/Linux/' PlotGenome; toc
% In Linux: tic; mcc -m -I './' -d './' PlotGenome; toc
% In Windows (32bit): tic; mcc -m -I 'D:\My Documents\Programming\MatLab\Work\' -d 'PlotGenome_eran'; toc
% In Windows (64bit): tic; mcc -m -I 'Z:\Tools\Matlab\Work\' -d 'Z:\IsoPlotter\' PlotGenome; toc
% 
% Note, this scirpt prints a complex object and uses the package export_figure for the printing.
% The package is avialable from: http://www.mathworks.com/matlabcentral/fileexchange/23629-exportfig
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Written by Eran Elhaik
% Written in : 11/13/12
% Ver : 2.00
% Website: http://code.google.com/p/isoplotter/
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ver 2.00 - cannot use chrlen, chr size is determined by the end of the chr.
%            Add workaround to make it print in linux.   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PlotGenome(input_file, output_file)

    disp('Start program PlotGenome');
    M1 = load(input_file);

    min_size = 0;
    iso_size = 300000;
    chr_num = numel(unique(M1(:,1)));

    chr_len_all = [];
    for chr=1:chr_num
        chr_len_all(chr) = max(M1(M1(:,1)==chr,3));
    end;

    %% Print chromsomal views

    for j=1:3
        genome_func = j; %1/2/3

        chr_pos_all = 1:-(1/(chr_num*2)):0;
        chr_len_all_prop = zeros(1,2*(chr_num));
        chr_len_all_prop(1:2:end) = chr_len_all./max(chr_len_all);

        chr = zeros(1,2*(chr_num));
        chr(1:2:end) = 1:chr_num;

        for chr_vec=1:numel(chr)

            curr_chr = chr(chr_vec);

            %Plot chromosome, not spaces
            if curr_chr>0
                chr_length = chr_len_all(curr_chr)/max(chr_len_all);
                disp([num2str(j) '. Now analyzing chr #' num2str(curr_chr)  ' with relative length of #' num2str(chr_length)]);

                Y_data = [chr_pos_all(chr_vec+1);chr_pos_all(chr_vec+1); chr_pos_all(chr_vec); chr_pos_all(chr_vec)];
                GenomePlot_func(genome_func, iso_size, min_size, M1(M1(:,1)==curr_chr,:), Y_data, chr_length);
                rectangle('Position',[0 chr_pos_all(chr_vec+1) chr_length chr_pos_all(chr_vec)-chr_pos_all(chr_vec+1)], 'LineWidth', 1, 'EdgeColor','k'); %draw a rectangle around the box to mark the edges
            end;
        end;
        fullscreen = get(0,'ScreenSize');
        axis([-0.001 1 0 1])
        grid off;
        box off;
        whitebg('w')
        axis off;
        set(gcf, 'Color', 'w');

        %Save figure
        filename = strfind(output_file, '.');
        if isempty(filename)
            disp('Output file must contain extension, e.g.: d:\file.tif');
            return;
        end;
        curr_output = [output_file(1:filename-1) '_' num2str(genome_func) output_file(filename:end)];
        disp(['Saving plot to ' char(curr_output) ' ...']);
        saveas(gcf, curr_output, 'tif');
%         print ('-dtiff', output_filename, res);
%         export_fig curr_output '-tif'
        close all;
    end;
    disp('End program PlotGenome');
end

