function [expVar, n, mu, p, F] = jANOVA1(data, groupIDs, dim, grps, GrandMeanMethod, CalcOmega2ExpVar)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [expVar, n, mu, p, F] = myANOVA1(data, groupIDs, dim, grps, GrandMeanMethod, CalcOmega2ExpVar)
% 
%   One-way ANOVA test/explained-variance analysis.  Calculates the proportion of total variance explained in 'data' array
%   by a given factor (which may have any number of levels or groups within in).  Also optionally returns F-statistics 
%   and associated p values.  (Note: this is essentially a stripped down, speeded up, and slightly enhanced version 
%   of Matlab's ANOVA1() -- the basic arguments work the same here as there).
%
%   Note (ADDED 4-15-14): Any data points with 0 data variance (SStotal = 0) now have expVar and F set = 0 (and p therefore = 1).
%     These are strictly undefined (divide by 0), and were previously left = NaN and then set = 0 post hoc.  
%     But I decided it makes most sense to just simply set these data points = 0 here, esp. since there might be other reasons for NaN values.
%
% Arg's:
%   data      Array of data to analyze (any arbitrary size -- analysis will be performed along dimension 'dim', and all other array dimensions
%              are treated as independent conditions [eg, different time points, trial windows, neurons, etc] to be analyzed separately)
%   groupIDs  Vector of group labels for each data point -- identifies which factor level/group each data point belongs to (eg, level 1,2,3, etc.)
%             Must be same length as size of DATA along dimension DIM to be analyzed, therefore groups must be same for all independent conditions
%             Currently, must be numeric, though could easily be coded to handle cell strings, for example.
%   dim       Dimension to perform ANOVA calc's along (all other dimensions treated as independent conditions to be analyzed separately)
%   grps      (optional): ID labels for each group expected to be in dataset 
%             (ie, should be output of unique(groupIDs), which is explicitly calc'd below if not input)
%   GrandMeanMethod   Method used to calculate grand mean for ANOVA formulas:
%                     0 = (default) Calc. grand mean as mean of all observations (arguably the more-standard ANOVA formula)
%                     1 = Calc. grand mean as mean of group means--causes less downward-biasing of expvars for unbalanced grp n's
%   CalcOmega2ExpVar  Formula for proportion of explained variance (see Snyder & Lawson 1993 Journal of Experimental Education or wikipedia entry for ANOVA):
%                     0 = (default) Standard eta-squared (R-squared) 'explained variance' = SSgrps / SStotal.  Biased upward, esp for small sample sizes.
%                     1 = Omega-squared. Bias-corrected explained variance = (SSgrps - dfgrps*MSerr)/(SStotal+MSerr) 
%
% Returns:
%   expVar    Proportion of total data variance explained by given factor, for each independent condition (non-'dim' dimension)
%             in data array, calculated as described above in 'CalcOmega2ExpVar' arg
%   n         # data points in each group
%   mu        Group means
%   p         p values corresponding to above expVar's, based on F-test
%   F         F statistics for above F-tests
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (nargin < 3) || isempty(dim)
  if isrow(data), dim = 2; 
  else            dim = 1; end
end
if (nargin < 5) || isempty(GrandMeanMethod),  GrandMeanMethod = 0;  end
if (nargin < 6) || isempty(CalcOmega2ExpVar), CalcOmega2ExpVar = false;  end

if (nargin < 4) || isempty(grps)              % If expected group labels were not input as arg, ...
  grps  = sort( unique(groupIDs(:)) );        % Find number of groups in dataset 
end                                           % (NOTE: currently does this for ALL data input, not separately along 'dim', so don't have diff grp's in diff slots there...)
nGrps = length(grps);

if length(groupIDs) ~= size(data,dim)
  error('GROUPIDS length (%d) must be same as size of DATA along dimension DIM (%d) to be analyzed', length(groupIDs), size(data,dim));
end

nDims           = ndims(data);
dataSize        = size(data);                 % Original size of data array
varSize         = dataSize;                   % Size of returned output variable arrays:           
varSize(dim)    = 1;                          %  'expVar', 'p', 'F' (reduced to singleton along analysis dim)

% Rearrange data array Y so that analysis (observations/trials) dim is dim 1, 
%  and everything else is concatenated into a vector (so arrays become size [#observations, total # indep. conditions])
if dim ~= 1
  dimPermute= [dim setxor([1:nDims],dim)];
  data      = permute(data, dimPermute);
end
if nDims > 2
  data      = reshape(data, dataSize(dim), []);
end
nIndepConds = size(data,2);                   % Number of independent conditions to analyze

% Calculate group means for each group (factor level) in dataset
[mu, n] = sub_calcGroupMeansAndNs(data, groupIDs, grps, nGrps, nIndepConds);

% If only 1 group in input data, or no data for any groups expected to be in data, give a warning and return NaNs for output vars
if (nGrps == 1) || any(n == 0),                 
  expVar= nan(varSize);                       
  if nargout > 3
    F   = nan(varSize);
    p   = nan(varSize);
  end
  % TODO: reshape mu here
  warning('Only 1 group in given data--returning NaNs'); 
  return;
end

% Calculate grand mean
N       = sum(n);                             % Total number of observations
if GrandMeanMethod == 0                       % Calculate grand mean as mean of all observations (standard ANOVA formula)
  grandMean = nanmean(data, 1);                   
else                                          % Calc. grand mean as mean of group means
  grandMean = mean(mu, 1);                    %  causes less downward-biasing of expvars for unbalanced grp n's
end

% Calculate explained variance for all data points
SSgrps    = sub_calcGrpsSumsOfSquares(mu,grandMean,n,nGrps,nIndepConds); % Groups Sum of Squares

SStotal   = nansum(abs(data - repmat(grandMean,[N 1])).^2, 1); % Total Sum of Squares (note: abs() is to deal w/ complex data)

if (nargout > 3) || CalcOmega2ExpVar
  dfGrps  = nGrps-1;                          % Groups degrees of freedom
  dfErr   = N-1 - dfGrps;                     % Error degrees of freedom
  MSerr   = (SStotal-SSgrps) ./ dfErr;      	% Error mean square  
end

if CalcOmega2ExpVar  
  expVar  = (SSgrps - dfGrps*MSerr) ./ ...  	% Omega-squared stat = bias-corrected explained variance cf. Snyder & Lawson, 1993; Olejnik & Algina, 2003; Buschman, et al, 2011
            (SStotal + MSerr);
  
else
  expVar  = SSgrps ./ SStotal;               	% Standard expvar -- Proportion of total data variance explained by factor groups ('eta-squared')
end
expVar(SStotal == 0) = 0;                     % expVar strictly undefined when no data variance (div by 0). Set = 0 for these cases (makes most sense) (ADDED 4-15-14)

% Calculate F-statistic and perform F-test to determine p value for all data points
if nargout > 3
  MSgrps  = SSgrps ./ dfGrps;                 % Groups mean square
  F       = MSgrps ./ MSerr;                  % F stat
  F(SStotal == 0) = 0;                        % F strictly undefined when no data variance (div by 0). Set = 0 for these cases (makes most sense) (ADDED 4-15-14)
  p       = 1 - fcdf(F,dfGrps,dfErr);         % P value for given F stat value
end


% Rearrange output variable arrays so in format/size expected based on original input arg's
if nDims > 2
  reshapeSize = [1 varSize(setxor(1:nDims,dim))];
  expVar  = reshape(expVar, reshapeSize);
  if nargout > 2
    mu    = reshape(mu, [nGrps reshapeSize(2:end)]);
  end
  if nargout > 3  
    F     = reshape(F, reshapeSize);
    p     = reshape(p, reshapeSize);
  end
end
if dim ~= 1
  expVar  = ipermute(expVar, dimPermute);
  if nargout > 2
    mu    = ipermute(mu, dimPermute);
  end
  if nargout > 3
    F     = ipermute(F, dimPermute);
    p     = ipermute(p, dimPermute);
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculate group means for each group (factor level) in dataset
function  [mu, n] = sub_calcGroupMeansAndNs(data, groupIDs, grps, nGrps, nIndepConds)  
  n             = nan(nGrps,1);
  mu            = nan(nGrps,nIndepConds);
  
  for iGrp = 1:nGrps                            % TODO: try to better 'Matlabize' this section?? (elim for loop)
    grpNdxs     = groupIDs == grps(iGrp);       % Find all observations in data belonging to given group
    n(iGrp)     = sum(grpNdxs);                 % #observations for given group
    mu(iGrp,:)  = nanmean(data(grpNdxs,:), 1);  % Group mean for given group
  end
  
  if nGrps == 1,                                % If only 1 group in input data, give a warning and return array
    mu = cat(1, mu, nan(1,nIndepConds));        %  padded with another 'group' of NaNs, as calling fun's may expect size of 2 grp's
  end  
  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculate groups Sums of Squares for all data points  
function  SSgrps = sub_calcGrpsSumsOfSquares(mu, grandMean, n, nGrps, nIndepConds)

  SSgrps    = zeros(1,nIndepConds);
  for iGrp = 1:nGrps                            
    SSgrps = SSgrps + n(iGrp).*abs(mu(iGrp,:) - grandMean).^2;  % Group Sum of Squares for given group (note: abs() is to deal w/ complex data)
  end   
  
  % Note: this code is actually ~an order of magnitude faster than the vectorized version below:
  % SSgrps2   = sum( repmat(n,[1 nInd]) .* (mu - repmat(grandMean,[nGrps 1])).^2, 1);
  