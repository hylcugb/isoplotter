%DjsSegmentation - DNA segmentation using Djs
% 
% Description: DjsSegmentation recursively partitions a sequence using Shannon
% entropy and a fixed threshold. 
% 
% Parameters: DjsSegmentation(dat, target_file, win_size, sizelim_, threshold)
%             dat - A txt file that contains the number of GC nucleotides
%                   per a certain window size (default is 32bp). DjsSegmentation accepts either the filename or a vector. 
%		    To create a file with the format, use the bash script: CalculateGCWindow.n.	
%             target_file -  Optional. Output filename. 
%             win_size - Optional. The window for which nucleotides were
%                        counted. Must be the same value as used to create the file. Default is 32bp.
%             size_lim_ - Optional. Minimum domain size. Should be a multiplication of win_size. Default is: 3008.  
%             threshold - a fixed value that halts the segmentation. Default is: 0.000058.  
% 
% Input: The program accepts both filename or a vector of GC content counts
%        per window size.
%           Example of input file of 32bp: 16 17 1 2 3 32 14 5 6
%           Example of input file of 1bp:  1  0  0 1 1 1  0  0 0 
% 
% Output file format: From base, To base, Domain length, Domain GC content,
%                     Domain GC content standard deviation
%  
% Examples: DjsSegmentation('Example.txt')
%           DjsSegmentation('Example.txt', 'output.txt')
%           DjsSegmentation('Example.txt', 'output.txt', 32, 3008)
%           DjsSegmentation('Example.txt', 'output.txt', 32, 3008, 0.0002)
%           DjsSegmentation([ 12 13 14 15 18 30 32 31 29 28])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Written by Eran Elhaik
% Written in : 07/01/08
% Ver : 1.20
% Website: http://code.google.com/p/isoplotter/
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function DjsSegmentation(dat, target_file, win_size, sizelim_, threshold)
	
if (nargin>0)
   if (nargin==0), disp('Not enough parameters! Format is: Source filename, Target filename (optional), window size (optional), minimum domain size (optional)'); return; end;
   if (nargin == 1), target_file = ['Djs_output.txt']; disp(['Output file is ' target_file]); win_size = 32; sizelim_ = 3008; threshold = 0.000058; end;
   if (nargin == 2), win_size = 32; sizelim_ = 3008; threshold = 0.000058; end;
   if (nargin == 3), sizelim_ = win_size^2*3; threshold = 0.000058; end;
   if (nargin == 4), threshold = 0.000058; end;
end

%If the dat file is a reference to a file
if ~isnumeric(dat)
    dat=load(dat);  %If this is a file link - load the file
end; 

%Create a row vector of data
dat=reshape(dat, [1 length(dat) ]);  

%Quality control of input sequence
if length([find(dat<0) find(dat>win_size)])
    disp(['Djs Error: Sequence should contain GC counts ranging from 0 to win_size. Exit program.']);
    return;
end;
    
%If sizelim is not a multiplication of win_size adjust it
sizelim = ceil(sizelim_/win_size);
if ~(sizelim==(sizelim_/win_size))
    disp(['Djs Warning: It is recommanded to set sizelim as a multiplication of win_size. System would adjust sizelim to ' num2str(sizelim*win_size)]);
end;

StopSegmentation = 0;
mat=[1; length(dat)];

%Start recursive segmentation of dat
while StopSegmentation<1,
    StopSegmentation=1;

    %Get fromto sizes (2 vecors with from and to)
    fromto = [mat(1:1:end-1) mat(2:1:end)];

    %Run on the length of the vecotr and for each border
    for i=1:size(fromto,1)
        parm = CalculateEntropy(dat(fromto(i,1):fromto(i,2))./win_size);

        %If the matrix is not empty, test if reached its limits
        if ((size(parm,2)) ~= 0)
            %If all conditions regarding threshold, min seg, and min sd
            if ( parm(2)>threshold && parm(1)>sizelim && (fromto(i,2)-parm(1)-fromto(i,1)>sizelim) )    
                mat = [mat; fromto(i,1)+parm(1)];
                StopSegmentation = 0;
            end;
        else
            disp(['Djs Error: Error on input sequence. ' num2str([i parm]) ]);
            disp(fromto(i,:));
        end;
    end;
    mat=sort(mat);
end;

%Finalize vector
fromto = [mat(1:1:end-1) mat(2:1:end)];         %Convert to 2 columns format
fromto32 = fromto*win_size; fromto32(:,1)=fromto32(:,1)+1; fromto32(1) = 1;    %Convert to Xbp format

%Calculate statistics for each domain: mean, std
iso_statistics = zeros(2,size(fromto32,1));
for stat = 1:size(fromto32,1)
    seq = dat(fromto(stat,1):fromto(stat,2))/win_size;  %sequence GC%
    iso_statistics(:,stat) = [mean(seq) std(seq)];   %This is std
end;
results = [fromto32 diff(fromto32,1,2)+1 iso_statistics'];

%Write results to file
fid = fopen(target_file, 'w+');
if fid
    fprintf(fid, '%d %d %d %.3f %.4f \n', results');
    fclose(fid);
else
    disp(['Djs Error: Cannnot open output file ' fout]);
end;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   CalculateEntropy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out=CalculateEntropy(gc);

%Segment length
len = length(gc);

%The sum of GC% content and AT% content, along the entire chr
sumgc = sum(gc); 
sumat = sum(1-gc);

%Mean GC and AT for the whole chr
FgcT = sumgc/len;
FatT = sumat/len;

% This return a single number. in here the sign is positive, it changes to
% negative after that. This line can potentially return NaN if one of the
% Fs is 0.
Htot = FgcT.*log2(FgcT) + FatT.*log2(FatT);

%Check the value of Htot
if ((Htot>0) | Htot<-1)
	display(['Djs Error: Serious Error with calculation of Htot. Check your sequence. ' num2str(Htot) ]);
    out = [];
    return;
end;

cx=cumsum(gc(1:len-1));
ax=cumsum(1-gc(1:len-1));
clear gc;

%Allocating memory
FgcS1 = cx./(1:len-1);                  %Cumulative sum divided by the segment index. In every index it is the GC% from 1-index.
FgcS2 = (sumgc-cx)./(len-1:-1:1);       %In every index it is the GC% from that index-end
FatS1 = ax./(1:len-1);
FatS2 = (sumat-ax)./(len-1:-1:1);
clear ax; clear cx;

%Make sure there will be no 0s for the log.
FgcS1(find(FgcS1==0))=0.001;
FgcS2(find(FatS1==0))=0.001;
FatS1(find(FgcS2==0))=0.001;
FatS2(find(FatS2==0))=0.001;

Hs1 = FgcS1.*log2(FgcS1) + FatS1.*log2(FatS1);
Hs2 = FgcS2.*log2(FgcS2) + FatS2.*log2(FatS2);
clear FgcS1; clear FgcS2; clear FatS1; clear FatS2;

%Check the value of Htot
if (sum(find(Hs1>0)) | sum(find(Hs1<-1)) | sum(find(Hs2>0)) | sum(find(Hs2<-1))   )
	display(['Djs Warning: Possible problem with calculation of Hs1 or Hs2 - correcting the problem']);
    Hs2(find(Hs2>0))=0;
    Hs2(find(Hs2<-1))=-1;
    Hs1(find(Hs1>0))=0;
    Hs1(find(Hs1<-1))=-1;
%	exit;
end;

%get fractures in the number of segmenrs from 0.1 to 1 and vice verse.
relen1 = (1:len)/len;     %this is the Li/L part
relen2 = (len-(1:len)+1)/len;

Djs = -Htot + relen1(1:end-1).*Hs1 + relen2(2:end).*Hs2;  %Final version

%Return the largest entropy value (Vmax) and its position (Dmax)
[Vmax Dmax] = max(Djs(1:end-1));

%Looks like for debugging purposes. the program gets 'i' all the time
out=[Dmax Vmax];

