%MapN1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       Print the Ns locations in a file
%
% This program accept 2 files, input & output. It reads the input and
% prints the Ns locations to a file. This script is faster than any perl
% script. The input files are the genome files. The script will produce
% several lines fo a Ns island that is in the size of several lines.
%
% Notice: The script assumes that the first line has been removed. If it
% wasn't than the whole mapping would be offset. The N mapping should be
% identical to the chr.apg description file. Sometimes is running out of
% memory.
%
% Notice: the windows may results in some of the sequence be lose. The N's are mapped precisely, so they may overlap the missing part.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Written by Eran Elhaik
% Written in : 11/12/05
% Ver : 17.00
% Website: http://code.google.com/p/isoplotter/
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ver 3.00: I corrected and optimized the code. it was incorrect before.
% Ver 4.00: I changed the printing to be precise
% Ver 5.00: I added a better documentation.
% Ver 6.00: I added a printing of seq without Ns, and 32bp. This version is
% different from previous one in its way of peration and memory allocation.
% This works and creates small djs files (with less spaces).
% Ver 7.00: Adding documentation. Simplifying the script and correcting bugs. I removed the division to 3, the new servers can handle that memory issue.
% I tested all possible N cases.
% Ver 8.00: Adding special case for header line. The script checks whether there is a header or not and ignores it if it does.
% Ver 9.00: Automatisizing for multiple files. Add QC for non ACTGN
% characters. Fixed a bug in the 32 calculation (last line is dropped) if
% not full.
% Ver 10.00: correct for filenames with multiple dots.
% Ver 11.00: I convert the chars to numbers to reduce memory size. 
% Input size up to 60Mb were tested on the computer. Tested well. but i had to remove the part that prints the sequence without N's.
% Ver 12.00: Add printing of original sequence size, and QC test so the
% final seq won't be larger than the original seq. Fixing bug in counting
% 32bp windows. Empty lines affected the count. Add a correction for small
% n's in the procedure that removes them. Adjusting script to deal with Ns
% erase
% Ver 13.00: add a variable to count the number of chars in the sequence
% line. Remove empty columns if there are any. I added print messages.
% Ver 14.00: Fixing a bug, if the last nucleotides are Ns and they are not included in the last window, then the addNs script would result in an error.
% Here, i trimmed the Ns domains depends on the window size
% Ver 15.00: When a non ACTGN character is found the program places N in there and notify the user.
% Ver 16.00: change target_filename to allow fasta files witout extension
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% input_dir = '../Genomes/bee/scaffolds_no_ns/'
% output_dir_seqs = '../IsoPlotter/bee/scaffolds_no_ns/'
% output_dir_32bp = '../IsoPlotter/bee/Djs/'
% MapN1('./Test/', './1/', './2/')

% input_dir=original_seq_dir;
% output_dir_seqs=List_ns;
% output_dir_32bp=Seq_bp_dir;

function MapN1(input_dir, output_dir_seqs, output_dir_32bp)

    program_name = 'MapN1';
    version = '15.00';
    tic;
    disp(['Started ' program_name ' version (' version ').']);
    win_size = 32;
    
    FilesList = dir(input_dir);
    FilesList(1:2) = [];

    disp(['#' num2str(size(FilesList,1)) ' files were found' ]);

    for files=1:size(FilesList,1)
%    for files=1:1

        %Get file name
        source_file = FilesList(files).name;
%        extension_size = length(source_file)-max(findstr(source_file,'.')); %was removed in 16.00
%        validation_ext = source_file(length(source_file)-extension_size : length(source_file)); %was removed in 16.00
%        target_filename = source_file(1:length(source_file)-(extension_size+1)); %was removed in 16.00
        if isempty(max(findstr('.',source_file))) 
            target_filename = source_file; 
        else
            target_filename = source_file(1:max(findstr('.',source_file))-1);
        end; %New in 16.00
        disp(['Analyzing: ' source_file ' File No.' num2str(files)]);

        seq_temp_temp = [];
        seq_temp_temp = textread([input_dir source_file],'%s','delimiter','\n','whitespace','','bufsize', 10095000);
        %disp(['Converting to Char . . . ']);
        seq_temp = char(seq_temp_temp);
        clear seq_temp_temp;
        disp(['    Finished Reading file: ',source_file ]); %Cannot tell what is the file size at this point because of the header line

        %Check if first line consits of a header
        if seq_temp(1)=='>'
            seq_temp(1,:)=[];
            line_size =  sum(seq_temp(1,:)~=' '); %sometimes there are empty spaces after the sequence which messes the size count
            seq_size = (size(seq_temp,1)-1) * line_size + sum(int8(seq_temp(end,:))~=win_size);
            disp(['    Header line was removed, sequence size is: #' num2str(seq_size) ]);
        end;

        % Remove empty columns
        if size(seq_temp,2)>line_size %If there are empty columns
            seq_temp(:,line_size+1:size(seq_temp,2))=[];
        end;
        
        %Check that all files fit ACTGN
        if sum(sum(~(seq_temp=='A' | seq_temp=='C' | seq_temp=='G' | seq_temp=='T'  | seq_temp=='N' | seq_temp=='a' | seq_temp=='c' | seq_temp=='g' | seq_temp=='t'  | seq_temp=='n' | seq_temp==' ')))
            disp('MapN1 Error: Found a weird character, replacing them with Ns')
            x=find (seq_temp~='A' & seq_temp~='C' & seq_temp~='G' & seq_temp~='T'  & seq_temp~='N' & ...
            seq_temp~='a' & seq_temp~='c' & seq_temp~='g' & seq_temp~='t'  & seq_temp~='n' & seq_temp~=' '); %Added in ver 15.00
            
            seq_temp(seq_temp~='A' & seq_temp~='C' & seq_temp~='G' & seq_temp~='T'  & seq_temp~='N' & ...
            seq_temp~='a' & seq_temp~='c' & seq_temp~='g' & seq_temp~='t'  & seq_temp~='n' & seq_temp~=' ')='N'; ; %Added in ver 15.00
            %continue;
        end;

	%calculate the final sequence size considering the part not covered by the window 
%	seq_size_no_ns = seq_size-sum(sum(upper(seq_temp)=='N')); %Final size, before trimming to windows
%	final_seq_size = floor(seq_size_no_ns/win_size)*win_size; %This is the final size after trimming to windows, 
	final_seq_size = floor(seq_size/win_size)*win_size;
        
        %Calculate GC%
        GC_num_line = (sum(upper(seq_temp)=='C') + sum(upper(seq_temp)=='G'));
        
        %Find Ns and calculate differences along the locations
        disp('    Finding N indexes...');
        find_ns_indexes =  find(upper(seq_temp)' == 'N');
        num_of_ns_before = length(find_ns_indexes); %How many Ns are initially

        find_ns_indexes = find_ns_indexes(find_ns_indexes<=final_seq_size); %Save only Ns in the size allowed to be saved; Ver 14.00
        
        num_of_ns = length(find_ns_indexes); %How many Ns are there
	num_of_ns_removed = num_of_ns_before-num_of_ns;

        %Save indexes of where Ns appear
        disp(['    Complete index findings: ' num2str(num_of_ns) ' Ns were found']);

        if isempty(find_ns_indexes)
            disp('    No Ns were found.');
        else
            %Saves start-end of NS appearnces
            Ns_mat = [0 find(diff(find_ns_indexes)>1)' length(find_ns_indexes)];

            %Calculate a 2 columns matrix in a fromto Ns format 
            %These are not the Ns borders, only their indexes in the find_ns_index vec
            fromto = [(Ns_mat(1:1:end-1)+1)' Ns_mat(2:1:end)']; %Notice, the size of each N island is: diff(fromto,[],2)+1 

            %Get the fromto in Ns indexes in the sequence
            Ns_indexes = find_ns_indexes(fromto);
            clear find_ns_indexes;
            
            %Print fromto indexes
            disp(['    Found ' num2str(sum(diff(Ns_indexes')+1)) ' Ns']);

            %check whether there is last incomplete line with Ns that should not be counted
%            if ~isinteger(size(seq_temp(:),1)/win_size)
%                %Remove Ns that are in the last incomplete line that will be removed anyway
%                disp('    Remove Ns from the last line...')
%                Ns_indexes(Ns_indexes>win_size*(floor(size(seq_temp(:),1)/win_size)-1)) = [];
%            end;
            
            if ~isempty(Ns_indexes)
                %%%%%%%%%%%%%%%%%%%%%%%%%
                %Print Ns coordinates
                %%%%%%%%%%%%%%%%%%%%%%%%%
                %Write results to file
                disp(['Writing Ns mapping to file ', target_filename]);
                fid1 = fopen([output_dir_seqs target_filename '.txt'], 'w+');
                fprintf(fid1, '%d %d \n', Ns_indexes');
                fclose(fid1);
                disp('  Clear memory');
                clear fromto;
                clear Ns_indexes;
                clear Ns_mat;
            end;
        end;

        %%%%%%%%%%%%%%%%%%%%%%%%%
        %Convert letters to numbers to save on memory
        %%%%%%%%%%%%%%%%%%%%%%%%%
        char_num = sum(sum(~isnan(seq_temp))); %overall size including empty ones
        
        seq_temp(upper(seq_temp)=='C')='1';
        seq_temp(upper(seq_temp)=='G')='1';
        seq_temp(upper(seq_temp)=='T')='0';
        seq_temp(upper(seq_temp)=='A')='0';
        seq_temp(upper(seq_temp)=='N')='2'; %Change N's to 2
        seq_temp(upper(seq_temp)==' ')='3'; %Change empty line to 3        
        seq_temp=int8(seq_temp)-48;
        
        %Check again the sequence
        if sum(sum(seq_temp>3 | seq_temp<0))
            disp('MapN1 Error: Found non ACGTN characters. Exiting.');
            return;
        end;
        
        %Check the number of Ns before and after the conversion
        if  sum(sum(seq_temp==2))~=(num_of_ns+num_of_ns_removed)
            disp('MapN1 Error: Unmatching number of Ns. Exiting.');
            return;
        end;
            
        
        %%%%%%%%%%%%%%%%%%%%%%%%%
        %Removing all N's
        %%%%%%%%%%%%%%%%%%%%%%%%%
        disp('    Removing all Ns and empty cells forming a single line');
        seq_temp = reshape(seq_temp', 1, char_num); %Make one line of seq
        seq_temp(seq_temp==2)=[]; %Removing all N's 
        seq_temp(seq_temp==3)=[]; %Removing all empty chars
        if(find (seq_temp == 2)>0), disp('MapN1 Error: ERROR! Ns were not removed!'); end;

        %Make sure GC% counts hold
        if sum([sum(GC_num_line)'- sum(seq_temp)'])
            disp('MapN1 Error: Error counting GC%. Exiting.');
            disp([sum(GC_num_line) sum(sum(seq_temp)) ]); 
        end;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%
        %Calculate win_size windows
        %%%%%%%%%%%%%%%%%%%%%%%%%
        % disp('Calculating win_size windows');
        
        %Check that after removing N's the sequence size remain the same
        if seq_size~=(length(seq_temp)+num_of_ns+num_of_ns_removed)
            disp('MapN1 Error: Error in cleaned sequence size. Exiting.');
            disp([seq_size length(seq_temp)+num_of_ns])
            return;
        end;

        %Chop last bp that don't fit in win_size and reshape lines to be of
        %win_size size
        disp(['    Size of seqs without Ns is: ' num2str(length(seq_temp))]);
        seq_temp = seq_temp(1:(floor(length(seq_temp)/win_size)*win_size)); %Losing the last bp that do not divide by win_size bp
        seq_temp = reshape(seq_temp, win_size, length(seq_temp)/win_size)'; 

	%Display the new sequence size if it is smaller than the final sequence
	seq_size_after_win = size(seq_temp,1)*size(seq_temp,2)+num_of_ns+num_of_ns_removed;
	if (seq_size-seq_size_after_win)
		disp(['Original sequence size # ' num2str(seq_size) ' size after windows #' num2str(seq_size_after_win) ' diff # ' num2str(seq_size-seq_size_after_win) ]);
	end;
	

        %Check that the overall size of the win_size windows is smaller than
        %the original sequence size
        if (seq_size-num_of_ns-num_of_ns_removed)-(size(seq_temp,1).*size(seq_temp,2))>win_size
            disp('MapN1 Error: win_size  windows>sequence size . Exiting.');
            disp([seq_size (size(seq_temp,1).*size(seq_temp,2))])
            return;
        end;
            
        %Save win_size bp sums in a single vec
        res = sum(seq_temp==1,2);

        %Check that the GC% is less than win_size  han the overall
        if sum(GC_num_line)-sum(res)>win_size
            disp('MapN1 Error: missing GC% in win_size bp windows GC%. Exiting.');
            disp([sum(GC_num_line) sum(res)])
            return;
        end;
        
        %Save in win_size bp form
        disp(['Writing win_size bp to file ', target_filename]);
        fid3 = fopen([output_dir_32bp target_filename '.txt'], 'w+');
        fprintf(fid3, '%d ', res);
        fclose(fid3);
    end; %end files
    disp(['End ' program_name ' in ' num2str(toc) ' minutes.']);
end
