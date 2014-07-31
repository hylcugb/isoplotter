% OrganizeIsoPlotter_scaffolds
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   Get segmentation results to a Single file (for scaffolds)
%
% The script reads segmentation files for all chr for a species, it then
% prints them all to a single file with the chr number in the first row.
% Based on OrganizeIsoPlotter.m 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Notice: it takes the script a few minutes to run.
% Website: http://code.google.com/p/isoplotter/
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function OrganizeIsoPlotter_scaffolds(input_dir, output_filename)
    program_name = 'OrganizeIsoPlotter_scaffolds';
    version = '1.00';
    tic;
    disp(['Started ' program_name ' version (' version ').']);

    %% Read all seqs from dir
    %Reading all files in the genome
    GenomeFilesList = dir(input_dir);
    GenomeFilesList(1:2) = [];

    Ma = [];
    Mb = [];
    scaffold_name = {};
    disp(['    Merging #' num2str(size(GenomeFilesList,1)) ' files ...']);
    for i = 1:size(GenomeFilesList,1)

        %Get file name
        curr_scaffold = GenomeFilesList(i).name;
        scaffold_name = curr_scaffold(1:strfind(curr_scaffold, '.')-1);
        
        M1 = load([input_dir curr_scaffold]);
        Ma = [Ma M1' ];

        Mb = [Mb (repmat({scaffold_name}, size(M1,1), 1))'];
    end;
    Ma = Ma';
    
    % Write results to file
    disp(['    Writing results in ' output_filename]);
    fid = fopen(output_filename, 'w+');
    for i=1:numel(Mb)
        fprintf(fid, '%s\t %-8d\t %-8d\t %-8d\t %.3f %.4f %d\n', char(Mb(i)), Ma(i,:));        
    end;
    fclose(fid);

    disp('End of OrganizeIsoPlotter_scaffolds program. ');
end
