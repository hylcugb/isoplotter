%RunIsoPlotterPipeline
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               Run Algo IsoPlotter
%                               --------------------
%
% Description : This script runs two scripts: MapN1 and IsoPlotter on all
% genomes. MapN1 creates Ns mapping, calculate 32bp Djs files, and creates
% chr files without any Ns (disabeled). After each chr is processed
% IsoPlotter segments the file
%
% Input : SourceDir direcotires
%
% Output: prints to dir segmentation results
%         
% Activate by : run it.
% Running time: ~5 minutes for chr
% Create exe file (2.5Min): tic; mcc -m -I 'Z:\Tools\Matlab\Work\' -d 'Z:\IsoPlotter\' RunIsoPlotterPipeline; toc
% Create exe file (2.5Min): tic; mcc -m -I 'D:\My Documents\Programming\MatLab\Work\' -d 'D:\Temp\' RunIsoPlotterPipeline; toc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Written by: Eran Elhaik
% Written date: 4/25/07
% Version : 5.00
% Website: http://code.google.com/p/isoplotter/
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ver 2.00: I added documentation and the call to MapN1 function
% Ver 3.00: Correct file filename with multiple dots
% Ver 4.00: Changes variable names, determined a conventional directory naming. Defined 3 variables.
% Ver 5.00: Remove internal folder reference, turn into a fully automated function. Add a
% call to OrganizeIsoPlotter.m. Add QC tests for input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% input_dir = 'Z:\nexsan\Genomes\Test\'; %for chrs
% output_dir = 'Z:\nexsan\Genomes\Test_Results\';
% input_dir = 'Z:\nexsan\Genomes\Test2\'; %for scaffolds
% output_dir = 'Z:\nexsan\Genomes\Test2_Results\';
% RunIsoPlotterPipeline Z:\Genomes\Human Z:\IsoPlotter\Human
%RunIsoPlotterPipeline Z:\nexsan\Genomes\Test\ Z:\nexsan\Genomes\Test_Results\)
function RunIsoPlotterPipeline(input_dir, output_dir, domain_min_size, ns_domain_min_size, win_size)
    warning('off', 'MATLAB:MKDIR:DirectoryExists');
    program_name = 'RunIsoPlotterPipeline';
    version = '5.00';
    tic;
    disp(['Started ' program_name ' version (' version ').']);

    if nargin==2
        domain_min_size = 3008; %default smallest compositional domain size
        ns_domain_min_size = 50000; %default ns island domain name
        win_size = 32; %default windows size to measure GC%
    elseif nargin==3 
        ns_domain_min_size = 50000; %default ns island domain name
        win_size = 32; %default windows size to measure GC%
    elseif nargin==4
        win_size = 32; %default windows size to measure GC%
    end;
    
    %QC folder size
    if ~exist(input_dir, 'dir')
        error('Input dir must be an existing directory.');
    end;
    if ~exist(output_dir, 'dir')
        error('Output dir must be an existing directory.');
    end;

    if input_dir(end) ~= filesep
        input_dir = [input_dir filesep];
    end;
    if output_dir(end) ~= filesep
        output_dir = [output_dir filesep];
    end;

    %Define directories
    output_List_ns_dir = (['1.List_ns' filesep]);
    output_Seq_bp_dir = (['2.Seq_32bp' filesep]);
    output_IsoPlotter_dir = (['3.IsoPlotter_no_ns' filesep]);
    output_H_test_dir = (['4.IsoPlotter_no_ns_H' filesep]);
    output_H_test_ns_dir = (['5.IsoPlotter_ns_H' filesep]);
    
    disp('Creating subfolders...');
    mkdir(output_dir, output_List_ns_dir);
    mkdir(output_dir, output_Seq_bp_dir);
    mkdir(output_dir, output_IsoPlotter_dir);
    mkdir(output_dir, output_H_test_dir);
    mkdir(output_dir, output_H_test_ns_dir);

    disp(['Calling MapN1 for N''s mapping and a preparation of GC proportion file...']);    
    MapN1(input_dir, [output_dir output_List_ns_dir], [output_dir output_Seq_bp_dir]);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Read all sequences from genome dir
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Reading all files in input dir
    FilesList = dir([output_dir output_Seq_bp_dir]);
    FilesList(1:2) = [];

    disp(['Partitioning #' num2str(size(FilesList,1)) ' sequences using IsoPlotter and testing for homogeneity... ']);
    %Going through chromosomes/files in the species directory
    for seq_range = 1:size(FilesList,1)
        %Get filenames and dir names
        source_file = FilesList(seq_range).name;
        ext_size = length(source_file)-max(findstr(source_file,'.'));
        source_file_fa = [source_file(1:length(source_file)-ext_size) 'fa'];
        IsoPlotter_res_dir = [output_dir output_IsoPlotter_dir source_file];

        %Perform analysis
        disp([num2str(seq_range) '.    Running IsoPlotter for ' [output_dir output_Seq_bp_dir source_file ]]);    
        IsoPlotterSegmentation([output_dir output_Seq_bp_dir source_file], [output_dir output_IsoPlotter_dir source_file], win_size, domain_min_size);

        disp([num2str(seq_range) '.    Running Homogeneity test for ' [output_dir output_Seq_bp_dir source_file ]]);    
        TestingHomogeneityStatistics([output_dir output_Seq_bp_dir source_file], [output_dir output_IsoPlotter_dir source_file], win_size, [output_dir output_H_test_dir source_file]);

        disp([num2str(seq_range) '.    Running Add N''s script for ' [output_dir output_Seq_bp_dir source_file ]]);    
        IsoPlotterAddNsH([output_dir output_H_test_dir], source_file, [output_dir output_List_ns_dir source_file], ns_domain_min_size, [output_dir output_H_test_ns_dir])
        disp([num2str(seq_range) '. Analyzing this sequence was completed in ' num2str(toc) '.']);
    end;
    disp(['Completed partitioning in ' num2str(toc) ' minutes.']);

    disp('Summarizing all results into a single file...')
    OrganizeIsoPlotter([output_dir output_H_test_dir], [output_dir  'IsoPlotter_no_ns_H.txt']);
    OrganizeIsoPlotter([output_dir output_H_test_ns_dir], [output_dir  'IsoPlotter_ns_H.txt']);
    
    disp(['End ' program_name ' in ' num2str(toc) ' minutes.']);
end
















