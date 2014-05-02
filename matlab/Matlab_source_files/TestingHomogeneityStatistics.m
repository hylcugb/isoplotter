%IsoPlotterHomogeneity - Calculating domain homogeneity in reference to the sequence
% Source: www
% 
% Description: IsoPlotterHomogeneity reads the original sequence and
% IsoPlotter segmentation results. It then determines whether each domain
% is homogeneous compared to the sequence.
%
% Parameters: IsoPlotterHomogeneity (input_file_seq, input_file_iso, win_size, output_file_seq)
%             input_file_seq - A txt file that contains the number of GC nucleotides
%                              per a certain window size (default is 32bp). IsoPlotterSegmentation accepts either the filename or a vector. 
%             input_file_iso - IsoPlotter segmentation results file.  
%             win_size - Optional. The window for which nucleotides were counted. Default is 32bp.
%             output_file_seq - Output Optional. Output filename. 
% 
% Input: the program accepts 2 files one that contains GC content counts
% per window sizw (1 2 3 6 32 0) and another that contains IsoPlotter
% segmentation results file [ 1 10 10 0.5 0.02; 11 20 10 0.3 0.01; ]
%
% Output file format: From base, To base, Domain length, Domain GC content,
%                     Domain GC content standard deviation, Homogeneous (1)
%                     or non homogeneous (0)
% 
% Examples: TestingHomogeneityStatistics('1.txt', 'H_1_result.txt')
%           TestingHomogeneityStatistics('1.txt', 'H_1_result.txt', 32)
%           TestingHomogeneityStatistics('1.txt', 'H_1_result.txt', 32, 'output.txt')
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Written by Eran Elhaik
% Written in : 7/12/08
% Ver : 2
% Website: http://nsm.uh.edu/~dgraur/eran/IsoPlotter/main.htm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ver 2.00 - correcting bug in case no significant domains are found
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function IsoPlotterHomogeneity(input_file_seq, input_file_iso, win_size, output_file_seq)
	
if (nargin>0)
   if (nargin<=1), disp('Not enough parameters! Format is: segmentation source filename, IsoPlotter results file, window size, Output file (optional)'); return; end;
   if (nargin == 2), win_size = 32; output_file_seq = ['H_ ' input_file_iso]; end;
   if (nargin == 3), output_file_seq = ['H_' input_file_iso]; end;
end

global F1;
  
sig=0.05;
segs_cutoff = 100;      %Notice, it is in REAL number. it will be divided later to get 32bp. 

non_h_segments = [];
h_segments = [];

%Load the chromosome data in Xbp
input_seq = load(input_file_seq);

%Quality control of input sequence
if length([find(input_seq<0) find(input_seq>win_size)])
    disp(['IsoPlotterHomogeneity Error: Sequence should contain GC counts ranging from 0 to win_size. Exit program.']);
    return;
end;

%Calculating parameters for the F-test
F1 = asin(sqrt(input_seq/win_size));	
ab = [length(F1)-1 var(F1)];
        
%Load IsoPlotter segmentation resutls
Djs_dat_comp_temp = load(input_file_iso);    

Djs_dat_comp = [floor(Djs_dat_comp_temp(:,1:3)./win_size) Djs_dat_comp_temp(:,4:end)]; %get Xbp borders, drop the last columns
Djs_dat_comp(1,1) = 1;
        
%Calculate domain lengths in Xbp
Djs_seq_length = floor((diff(Djs_dat_comp_temp(:,1:2),1,2)+1)./win_size); 

%Get segments >threshold
Filter_segments = floor(segs_cutoff/win_size);
Djs_long_domain_ind = find(Djs_seq_length > Filter_segments);
Domains_in_analysis = length(Djs_long_domain_ind)/length(Djs_seq_length);   %Single number - proportion of domains included in analysis

if (Djs_dat_comp(end,2) > length(input_seq)*win_size)
    disp('ERROR! Segmentation file size > chr file size');
    return;
end;

% Homogeneity Analysis - Detect nonhomogeneous domains
D_P = zeros(1,length(Djs_long_domain_ind));        

Segment_borders = Djs_dat_comp(Djs_long_domain_ind,:);  %Get real borders for the long domains
for i=1:length(Djs_long_domain_ind)
    [H P] = vartest_h(Segment_borders(i,:), ab, sig, 'left');
    D_P(i) = P;
end;
        
%Calculate FDR correction
m = length(D_P);              
alpha = sig;
[pi I] = sort(D_P,'ascend');
y = ((1:m)./m)*alpha;         
sig_domains = max(find(pi<=y));         

% Ver 2.00 - correction
%Reject all null hypotheses H(1) … H(k). so they are all H1 - more homogeneous than the chr
if ~isempty(sig_domains)
    disp('sig_domains is empty')
    Homogeneous_domains_in_analysis = find(I<=sig_domains)'; %These are the indexes of the long H domains --corrected--
    putative_h_domains = length(Homogeneous_domains_in_analysis);

    %Save indexes of nonhomogeneous domains
    under_domains_indexes = Djs_long_domain_ind(find(I>sig_domains)'); %everything larger than the sig is non-H
else
    Homogeneous_domains_in_analysis = [];
    putative_h_domains = 0;
    under_domains_indexes = numel(I);
end;

       
%Update segmentation results + H test results (before concatenating)
Djs_dat_comp_temp(:,6) = 0;
Djs_dat_comp_temp(Homogeneous_domains_in_analysis,6) = 1;

non_h_segments = length(under_domains_indexes)/length(Djs_long_domain_ind);
h_segments = putative_h_domains/length(Djs_long_domain_ind);

%     Check that all proportions sum to 1
if abs(non_h_segments+h_segments-1)>0.001 
    disp('Warning! Statistics do not sum to 1. Your sequence may be very short to qualify fot the test');
    disp([non_h_segments+h_segments])
end;

%Creating updated segmentation file with homogeneity test results
fid=fopen(output_file_seq,'w+');
fprintf(fid, '%-8d %-8d %-8d %.3f %.4f %d \n', Djs_dat_comp_temp');
fclose(fid);
        

% -----------------------------------
% Calculate F-test
function [h,p,ci,stats] = vartest_h(x,y,alpha,tail)

global F1;
df2=y(1);
vary=y(2);

x=F1(x(1):x(2));
y=F1;

dim = find(size(x) ~= 1, 1);
if isempty(dim), dim = 1; end

% If we haven't been given an explicit dimension, and we have two
% vectors, then make y the same orientation as x.
if isvector(x) && isvector(y)
	if (size(x,1)==1)
	    y = y(:)'; 
	else
	    y = y(:);
	end
end

if ~isscalar(dim) || ~ismember(dim,1:ndims(x))
    error('stats:vartest2:BadDim', ...
      'DIM must be an integer between 1 and %d.',ndims(x));
end        

% Make sure all of x's and y's non-working dimensions are identical.
sizex = size(x); sizex(dim) = 1;
sizey = size(y); sizey(dim) = 1;
if ~isequal(sizex,sizey)
    error('stats:vartest2:InputSizeMismatch',...
          'The data in a 2-sample F test must be commensurate.');
end

% Compute statistics for each sample
[df1,varx] = getstats(x,dim);
%[df2,vary] = getstats(y,dim);

% Compute F statistic
F = NaN(size(varx),superiorfloat(varx,vary));
t1 = (vary>0);
F(t1) = varx(t1) ./ vary(t1);
t2 = (varx>0) & ~t1;
F(t2) = Inf;

% Compute the correct p-value for the test, and confidence intervals
% if requested.
p = fcdf(F,df1,df2);
if nargout > 2
ci = cat(dim, repmat(0,size(F)), ...
	      F./finv(alpha,df1,df2));
end

% Determine if the actual significance exceeds the desired significance
h = cast(p <= alpha, class(p));
h(isnan(p)) = NaN; % p==NaN => neither <= alpha nor > alpha

if nargout > 3
    stats = struct('fstat', F, 'df1', cast(df1,class(F)), ...
                               'df2', cast(df2,class(F)));
    if isscalar(df1) && ~isscalar(F)
        stats.df1 = repmat(stats.df1,size(F));
    end
    if isscalar(df2) && ~isscalar(F)
        stats.df2 = repmat(stats.df2,size(F));
    end
end

% -----------------------------------
function [df,varx] = getstats(x,dim)
%GETSTATS Compute statisics for one sample
  
% Get sample sizes and df
xnans = isnan(x);
nx = sum(~xnans,dim);
df = max(nx-1, 0);

% Get means
x(xnans) = 0;
xmean = sum(x,dim) ./ max(1,nx);

% Get variances
if isscalar(xmean)
   xcntr = x - xmean;
else
   rep = ones(1,ndims(x));
   rep(dim) = size(x,dim);
   xcntr = x - repmat(xmean,rep);
end
xcntr(xnans) = 0;
varx = sum(abs(xcntr).^2,dim);
t = (df>0);
varx(t) = varx(t) ./ df(t);
varx(~t) = NaN;

% Make df a scalar if possible, better for later finv call
if numel(df)>1 && all(df(:)==df(1))
   df = df(1);
end

