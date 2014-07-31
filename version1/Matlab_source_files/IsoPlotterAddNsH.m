% IsoPlotterAddNsH
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    Update segmentation domain with Ns coordinates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This sciprt reads a segmentation results file and a Ns mapping file. It
% updates the indexes of the segmenation file according to the Ns and
% prints a new file.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Written by Eran Elhaik
% Written in : 7/2/08
% Ver : 10.00
% Website: http://code.google.com/p/isoplotter/
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ver 2.00: This version is a function that accepts parameters are writes
% the results to file.
% Ver 4.00: I corrected a serious bug in the mapping calculation resutled
% in negative sized domains. I added more controls and added return values
% to the function that would return an error. I corrected all bugs. now the bias is minimal
% Ver 5.00: After adding 1 to n mapping i had to adjust the a future border calculation.
% Ver 6.00: This version includes a cutoff value that can have the Ns
% include in the domain.
% Ver 7.00: Adjusted to win xp os. not working. Documenting. Add a
% coreection to bug - Ns as a last line were not included because sometimes
% there was a gap in the sequence. I added a_addings variable.Also, no need
% for two functions, if the size of Ns island included in the domain can be
% a variable it is possible to use this function as a general/specific
% case. I added 2 more parameters to the function
% Ver 8.00: Works without bugs, get the numbers exactly as in the file. the
% only bugs shown are for NCBI's N files that i did on Atlnatis because
% tehre were to big for MapN1 to work with. when i work with MapN2 that
% reads the NCBI files, there are some single Ns that are not covered and
% they cause the error.
% Ver 9.00: Optimizing according to Matlab. also still have some bugs. this
% version matches the new IsoPlotter format (1-2, 3-4)
% Ver 10.00: If no N file is found - output the input fiel and exit
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%seg_dir = 'D:\My Documents\University\Dan Lab\01. Isochores Project\Data\IsoPlotter\Pan_troglodytes\';
%seg_file = '1.txt_res';
%n_file = 'D:\My Documents\University\Dan Lab\01. Isochores Project\Data\IsoPlotter\Pan_troglodytes\Ns\1.fa';

% seg_dir = ('D:\Eran\IsoPlotter\Pan_troglodytes\');
% seg_file = ('11');
% n_file = ('D:\Eran\Djs\Pan_troglodytes\Ns\11.fa');
% n_cutoff = 0;
% target_dir = 'D:\Eran\'

% species = 'Homo_sapiens'
% seg_dir = ['../IsoPlotter/' species '/4.IsoPlotter_no_ns_H/'];
% seg_file = 'hs_alt_HuRef_chr15.txt';
% n_file = ['../IsoPlotter/' species '/1.List_ns/hs_alt_HuRef_chr15.txt'];
% n_cutoff = 0;
% target_dir = ['../IsoPlotter/' species '/5.IsoPlotter_ns_H/'];

function IsoPlotterAddNsH(seg_dir, seg_file, n_file, n_cutoff, target_dir)

    %disp(['IsoPlotterAddNs (' num2str(n_cutoff) ') Load files']);

    %Load segmentation and result files
%     disp('    Loading chr file...');
    a = load([seg_dir seg_file]);
    
    %Define target file.
    fout = ([target_dir '/' seg_file ]);

    %If there is no N file, output the input file and exit	
    if ~(exist(n_file, 'file'))
    	disp('   No N files was found, outputting the input file... ')

        %Write results to file
        fid = fopen(fout, 'w');
        fprintf(fid, '%-8d %-8d %-8d %.3f %.4f %d \n', a');
        fclose(fid);
    	
    	return;
    end;
    
%     disp('    N list file...');
    b = load(n_file);   %Load N file
%     b = b(find((b(:,2)-b(:,1)) > 0),:);	
    
    c = [];
    n_shift = 0;
    passed_ns = 0;
    mapping_index = 0;	
    mapping_ind = [];	
%     n_cutoff = 10000; %Ns of smaller sizes than that will be included in the domain
    a_addings = 0;
    
    %Run over the domain list and check for every domain
%    disp('IsoPlotterAddNs: Go over domain list');
    for i=1:size(a,1)
%     for i=1:304

        %Find Ns that their start index is within the current a(i,:) domain
        mapping = find( (b(:,1)>=(a(i,1)+n_shift)) & (b(:,1)<=(a(i,2)+n_shift)) );

        if (size(mapping,2)==1)
            mapping = mapping';
        end;
        
        if ~isempty(mapping)

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Find all Ns that needs to be included in calculation
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            %seq_length of a domain
            seq_length = a(i,2)-a(i,1)+1;

            %length of all Ns that were found so far in mapping
            all_ns = sum(diff(b(mapping,1:2),1,2)+1);

            %This is the location of the farthest border, including all Ns and
            %the domain length
            next_border = a(i,1)+n_shift + length(mapping) + seq_length + all_ns; %Summing all up

            %Check if we passed an N islad on the way - find Ns from the
            %location of the last N to the last border.
            temp = find(  (b(:,1)<=next_border) & (b(:,1)>b(max(mapping),2)));  %Correcting bug >=

            %If found - add them to mapping and repeat calcualtion
            while ~isempty(temp)
                if (size(temp,2)==1), temp = temp'; end ;
                mapping = [mapping temp];

                %length of all Ns that were found in mapping
                all_ns = sum(diff(b(mapping,1:2),1,2)+1);

                next_border = a(i,1)+n_shift + length(mapping) + seq_length + all_ns; %Summing all up

                %Check if we passed an N islad on the way
                temp = find( (b(:,1)<=next_border) & (b(:,1)> b(max(mapping),2)));
            end; %end while

            
%             fprintf(1,'\n');
%             disp([num2str(i) '. Found original Mapping: ' num2str(mapping) ]);
%             disp(['Current Ns : ']);
% %             disp([b(mapping,:)]);
%             disp(['Including all Ns is resized domain: '	 num2str([a(i,1)+n_shift next_border]) ]);

            mapping_index = max([mapping_index max(mapping)]);
            mapping = setdiff(mapping, mapping_ind);    %In case there are multiple mapping values
            mapping_ind = [mapping_ind mapping];

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            initial_a_index = a(i,1) + n_shift; %The starting domain index
            inc_counter = 0;    %Counter of Ns within the domain

            %Run on the domain with all Ns island found    
%             disp([num2str(i) '. Now adding ' num2str(size(mapping,2)) ' Ns']);
            for m=1:size(mapping,2)

                if (diff(b(mapping(m),:),1,2)+1 > n_cutoff)
    
%                     disp('Found Long Ns island');
                    %Calculate Pre N domains ([1 10] -> [1 5 GC%])
                    
                    %Check if b is not on start of a domain
                    if b(mapping(m),1)>initial_a_index+inc_counter
                        c = [c [initial_a_index+inc_counter b(mapping(m),1)-1 (b(mapping(m))-(initial_a_index+inc_counter)) a(i,4) a(i,5) a(i,6)]' ]; 
                    else
%                         disp('Special case, b opens the domain, skipping a');
                    end;

%     disp([num2str(i) '. Current a Domain: '	num2str(initial_a_index+inc_counter) ' ' num2str(b(mapping(m),1)) ]);

                    %Calculate N domains ([5 10 0%])
                    c = [c [b(mapping(m),1) b(mapping(m),2) (b(mapping(m),2)-b(mapping(m),1))+1 0 0 0]' ]; 

%     disp([num2str(i) '. Current N Domain index=' num2str(mapping(m)) ', '	num2str(b(mapping(m),1))  ' ' num2str(b(mapping(m),2)) ]);

                    %Calculate the shift of N's to position the next start-a border.
                    inc_counter = b(mapping(m),2) - initial_a_index+1;    %Distance over (a) from initial_a_index

                else
%                     disp(['b smaller than n_cutoff, skipping Ns ' ]);                    
                end; %End if
                
                
                %How many Ns have been used (count all) - save in a vec
                passed_ns = [passed_ns passed_ns(end)+1+(b(mapping(m),2)-b(mapping(m),1))];

                n_shift = passed_ns(end); %Total number of Ns that have passed

                %Test that all counting is correct
                if (n_shift ~= sum(diff(b(1:mapping(m),:),1,2)+1))
%                    disp(c');
%                    disp([i m n_shift sum(b(mapping(m),2)-b(mapping(m),1)) ]);
%                    disp([i m n_shift sum(b(1:mapping(m),2)+1-b(1:mapping(m),1)) ]);            

                    disp(' ***     ERROR! N counter is bad. Exit running   (1) . Pausing   ***');
%                     return;
                    pause
                end;
        
            end; %end for size(mapping,2)


            %Set first border from the last one
            if size(c)>0, last_border=c(2,end)+1; else last_border=1; end;
            
            %Special case that the b was all the way to the border
            if (last_border <= a(i,2)+n_shift)
                c = [c [last_border a(i,2)+n_shift (a(i,2)+n_shift-last_border)+1 a(i,4) a(i,5) a(i,6)]'];
            end;
            
% disp(['Closing Domain: '	num2str(last_border)  ' ' num2str(a(i,2)+n_shift) ]);

        %No Ns were found within a domain - shift domain indexes
        else
%             disp('Second close')
            c = [c [a(i,1)+n_shift a(i,2)+n_shift a(i,2)-a(i,1)+1 a(i,4) a(i,5) a(i,6)]' ];
        end; %end if

        
        if ( ~isempty(find(diff(mapping_ind)>1)) && (length(mapping_ind)>1))
            disp(' ***     Error! Skipped Mapping      ***');
            disp(mapping);
            return;
            pause
        end;
    
        %Check for negative domain length
        if (length(find(c(3,:)<0)))
            disp(' ***     ERROR! Negative Segments size!       ***');
            return;
        end;
        
        %Announce completion
        if (~mod(i,10000))
            disp(['Completed ' num2str(i) ' sequences - out of ' num2str(size(a,1))]);
        end;

        %Test countings
        if (n_shift ~= sum(diff(b(1:mapping_index,:),1,2)+1))
            disp(c');
            disp([i n_shift sum(b(mapping_index,2)-b(mapping_index,1)) ]);
            disp([i n_shift sum(b(1:mapping_index,2)+1-b(1:mapping_index,1)) ]);            
            
            disp(' ***     ERROR! N counter is bad. Exit running   (2)    ***');
            return;
        end;

        
        %Check that borders have consequecings numbers
        if (find(c(1,2:end)'-c(2,1:end-1)'>1))
            disp(' ***     ERROR! Borders do not follow - there are jumps ***');
            disp([num2str(length(find(c(1,2:end)-c(2,1:end-1)>1))) ' Borders jump! Terminate run. i=' num2str(i) ' Violation was found in lines: ']);
            disp(find(c(1,2:end)-c(2,1:end-1)>1));
            disp(c');
            disp(i);
%             return;
%             i
            pause
        end;    
% c'
% i
% pause

    end; %End main for loop [1:size(a,1)]
    
    
    %% If there are Ns island remaining
    if max(mapping_ind)<size(b,1)
        disp([max(mapping_ind) size(b,1)]);
        if (size(b,1)-max(mapping_ind) > 1)
            if (i<size(a,1)), disp('**** Check the i ****'); end;
            disp ('*****       More than one Ns island remaining!           ****'); 
            disp([size(b,1) max(mapping_ind) size(b,1)-max(mapping_ind)]);
%             return; 
        end;
        
        %If the distance between last member of c and last member of b is
        %larger than 1, then add an additional line so there won't be gaps
        if (b(size(b,1),1)-c(2,end)' > 1) 
            a_addings = b(size(b,1),1)-c(2,end)'-1;
            c = [c [c(2,end)' b(size(b,1),1) 1+b(size(b,1),1)-c(2,end)' c(4,end)' c(5,end)' c(6,end)']' ];
        end;
            c = [c [b(size(b,1),1) b(size(b,1),2) 1+b(size(b,1),2)-b(size(b,1),1) 0 0 0]' ];        
            
        mapping = size(b,1);
        mapping_ind = [mapping_ind mapping];
        passed_ns = [passed_ns passed_ns(end)+1+(b(mapping,2)-b(mapping,1))];
        n_shift = passed_ns(end); %Total number of Ns that have passed
    end;
    
    
    %Check that all is were used and if not display error
    if (i<size(a,1))
        disp('ERROR! Terminate early');
        fprintf(1,'\n');
        disp(['i= ' num2str(i) ' Current a domain: '	num2str(a(i,1)+n_shift) ' ' num2str(a(i,2)+n_shift) ]);
        return;
    end;
    
%     disp('IsoPlotterAddNs: Complete adding Ns');
    c = c';

    error = 0;
    %Control - check that the length of the new matrix is the old+Ns
    total_nucleotides = c(end,2)-a(i,2)-a_addings-n_shift;
    if (total_nucleotides  ~= 0 )                        %corrected with a_addings
        
        if total_nucleotides>10
            disp(' ***     ERROR! Final length is in error!     ***');
            disp([(c(end,2)-a(i,2)-a_addings) n_shift sum(diff(b,1,2)+1)]); 
        else
            disp(['--- Notice total_nucleotides is short in ' num2str(total_nucleotides) ' ----']);
        end;
    end;
    if (length(find(c(:,3)<0)))
        disp(' ***     ERROR! Negative Segments size!       ***');
        disp([i])       
        error = 1;        
    end;
    if (length(find(c(:,3)==0)))
        disp(' ***     ERROR! 0 size segments were found!       ***');
        disp([i]);
        error = 1;        
    end;
    if (sum(diff(([b(:,1) (b(:,2)+1) ]),1,2)) ~= n_shift)
    	disp(' ***     ERROR! Not all Ns were used in the sequence      ***');
        [sum(diff(([b(:,1) (b(:,2)+1) ])')) n_shift  sum(diff(([b(:,1) (b(:,2)+1) ])'))-n_shift ]
    	%Becuase of the way i created the borders there might be a small shift
        error = 1;        
    end;
    if ( (n_shift + a(end,2)) - c(end,2) > 0)
    	disp(' ***     ERROR! Adding Ns numbers dont match      ***')
        error = 1;        
    end;
    if (find(c(2:end,1)-c(1:end-1,2)>1))
        disp(' ***     ERROR! Borders to not follow - there are jumps ***');
        disp([num2str(length(find(c(2:end,1)-c(1:end-1,2)>1))) ' Borders jump']);
        [find(c(2:end,1)-c(1:end-1,2)>1)]
        error = 1;
    end;
    if ~error
%        disp('     IsoPlotterAddNs: Completed Successfuly');
    end;
  
    
  
	% Print b and its border length
	%     [b(:,1) b(:,2) b(:,2)-b(:,1)]

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Wait to check if everything is ok
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%Write results to file
	fid = fopen(fout, 'w');
	fprintf(fid, '%-8d %-8d %-8d %.3f %.4f %d \n', c');
	fclose(fid);
end
