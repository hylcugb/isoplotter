% OrganizeIsoPlotter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   Get segmentation results to a Single file
%
%The script reads segmentation files for all chr for a species, it then
%prints them all to a single file with the chr number in the first row.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Notice: it takes the script a few minutes to run.
% Website: http://code.google.com/p/isoplotter/
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ver 2.00: Adjusted for single species with multiple files that may or may not ave chr number.
% Ver 3.00: Adjusted to new directory structure of IsoPlotter results. Producing 2 results file with and without Ns.
% Ver 4.00: Adjust to different filenames that are either chromosomal or scaffolds
% Ver 5.00: Convert to a fucntion, integrate into the pipeline
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function OrganizeIsoPlotter(input_dir, output_filename)
    program_name = 'OrganizeIsoPlotter';
    version = '5.00';
    tic;
    disp(['Started ' program_name ' version (' version ').']);

    %% Read all seqs from dir
    %Reading all files in the genome
    GenomeFilesList = dir(input_dir);
    GenomeFilesList(1:2) = [];

    Ma = [];
    disp(['    Merging #' num2str(size(GenomeFilesList,1)) ' files ...']);
    for genomes=1:size(GenomeFilesList,1)

        %Get file name
        genome_source_file = GenomeFilesList(genomes).name;

        %If this is a chromosomal file, it would contain the string 'chr':
        if findstr(genome_source_file,'chr')
            chr_name = genome_source_file(findstr(genome_source_file,'chr')+3:max(findstr(genome_source_file,'.'))-1);
            chr_num = str2num(chr_name(min(regexp(chr_name, '[0-9]')):max(regexp(chr_name, '[0-9]')))); %convert to numeric
            disp(['        Merging : ' genome_source_file ' file #' num2str(chr_num) '.']);

        %if this is a scaffold, simply number the files according to their order in the directory
        else
            chr_num = genomes;
        end; 

        M1 = load([input_dir genome_source_file]);
        Ma = [Ma ([ repmat(chr_num, size(M1,1), 1) M1])' ];
        clear M1;
    end;

    Ma = sortrows(Ma',[1 2])';

    %Write results to file
    disp(['    Writing results in ' output_filename]);
    fid = fopen(output_filename, 'w+');
    fprintf(fid, '%-1d\t %-8d\t %-8d\t %-8d\t %.3f %.4f %d\n', Ma);        
    fclose(fid);

    disp('End of OrganizeIsoPlotter program. ');
end
