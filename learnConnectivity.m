function connWeights = learnConnectivity(varargin)
%learnConnectivity Learns connectivity maps over one or two regions
%   learnConnectivity estimates function connectivity maps from fMRI data,
%   as described in these two papers:
%
%   C. Baldassano, M.C. Iordan, D.M. Beck, L. Fei-Fei. "Voxel-Level
%   Functional Connectivity using Spatial Regularization." Neuroimage.
%   2012 Jul 28;63(3):1099-1106. doi: 10.1016/j.neuroimage.2012.07.046
%
%   C. Baldassano, M.C. Iordan, D.M. Beck, L. Fei-Fei. "Discovering
%   Voxel-Level Functional Connectivity Between Cortical Regions."
%   Machine Learning and Interpretation in Neuroimaging Workshop, Neural
%   Information Processing Systems (NIPS) 2012.
%
%   Please cite these papers if you publish results using this function.
%   Any comments, questions, or bug reports can be directed to
%   Chris Baldassano <chrisb33@cs.stanford.edu> <chrisbaldassano.com>
%
%   Using learnConnectivity in 'both' mode (described below) requires 
%   the CVX package to be installed. CVX can be downloaded for free from
%   http://cvxr.com/cvx/ (this function has been tested with CVX v1.22 and 
%   v2.0 beta).
%
%   This function has been tested on MATLAB versions 2010b, 2011a, 2011b,
%   and 2012a.
%
%
%
%   All function arguments are specified as property, value pairs. The
%   function has two main modes of operation:
%
%   learnConnectivity('type','one', ...)
%   In this mode, a connectivity map is learned over region 1 using the
%   mean timecourse from region 2 as the seed. The required inputs are:
%   adj1: adjacency matrix for region 1 (numvox1 x numvox1)
%   bold1: timecourses for region 1 (numvox1 x numTimepoints)
%   bold2: timecourses for region 2 (numvox2 x numTimepoints)
%   The smoothness parameter lambda can be selected by either:
%   a) passing in a value for lambda
%   b) automatically selecting a value for lambda using cross-validation.
%      In this case, a vector runLabels (1 x numTimepoints) must be
%      specified. During cross-validation, the test set will consist of all
%      timepoints labeled 1, then all those labeled 2, then all those
%      labeled 3, etc.
%   Optional parameters are:
%   allowNegative: a 0/1 flag specifying whether negative weights are
%                  allowed. If 0, then weights are constrained to be
%                  positive (using quadprog rather than a closed-form
%                  equation) (default: 1)
%   lambdaVector: a row vector of all lambda values to try during
%                 cross-validation (default: logspace(-6,8,60))
%   zscore: a 0/1 flag specifying whether voxels should be normalized to
%           have standard deviation of 1 (default: 0)
%   quiet: a 0/1 flag specifying whether warnings and other output should
%          be suppressed (default: 0)
%
%   The return value is a column vector of length numvox1 giving the
%   relative connectivity strength for each voxel (weights may be positive
%   or negative). This connectivity map is computed using the closed-form
%   solution to the quadratic objective rather than CVX.
%
%
%
%
%   learnConnectivity('type','both', ...)
%   In this mode, connectivity maps are learned simultaneously over regions
%   1 and 2. The required inputs are:
%   adj1: adjacency matrix for region 1 (numvox1 x numvox1)
%   adj2: adjacency matrix for region 2 (numvox2 x numvox2)
%   bold1: timecourses for region 1 (numvox1 x numTimepoints)
%   bold2: timecourses for region 2 (numvox2 x numTimepoints)
%   lambda: smoothness parameter
%   Optional parameters are:
%   zscore: a 0/1 flag specifying whether voxels should be normalized to
%           have standard deviation of 1 (default: 0)
%   quiet: a 0/1 flag specifying whether warnings and other output should
%          be suppressed (default: 0)
%   numRestarts: number of initializations to try (default: 20)
%   trustTheta: constraint tolerance during optimization (default: 0.05)
%   trustDelta: trust radius during optimization (default: sqrt(0.05))
%   stepEpsilon: convergence criteria (default: sqrt(0.05)/100)
%   matchEpsilon: criteria for matching redundant solutions (default: 0.01)
%
%   The return value consists of two cells corresponding to regions 1 and
%   2. Each cell contains a matrix of size numvox x numSolutions, i.e.
%   connWeights{1}(:,i) and connWeights{2}(:,i) give the connectivity maps
%   over regions 1 and 2 for solution i. The connectivity weights are all
%   nonnegative.


% Copyright (c) 2015, Christopher Baldassano, Stanford University
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
%     * Redistributions of source code must retain the above copyright
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright
%       notice, this list of conditions and the following disclaimer in the
%       documentation and/or other materials provided with the distribution.
%     * Neither the name of Stanford University nor the
%       names of its contributors may be used to endorse or promote products
%       derived from this software without specific prior written permission.
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
% ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
% WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
% DISCLAIMED. IN NO EVENT SHALL CHRISTOPHER BALDASSANO BE LIABLE FOR ANY
% DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
% (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
% LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
% ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
% (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
% SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


inpParse = inputParser;
inpParse.FunctionName = 'learnConnectivity';
inpParse.addParamValue('type',[]);
inpParse.addParamValue('adj1',[]);
inpParse.addParamValue('adj2',[]);
inpParse.addParamValue('bold1',[]);
inpParse.addParamValue('bold2',[]);
inpParse.addParamValue('runLabels',[]);
inpParse.addParamValue('lambda',0);
inpParse.addParamValue('lambdaVector',logspace(-6,8,60));
inpParse.addParamValue('numRestarts',20);
inpParse.addParamValue('zscore',0);
inpParse.addParamValue('trustTheta',0.05);
inpParse.addParamValue('trustDelta',sqrt(0.05));
inpParse.addParamValue('stepEpsilon',sqrt(0.05)/100);
inpParse.addParamValue('quiet',0);
inpParse.addParamValue('matchEpsilon',0.01);
inpParse.addParamValue('allowNegative',1);

inpParse.parse(varargin{:});
type = inpParse.Results.type;
adj1 = inpParse.Results.adj1;
adj2 = inpParse.Results.adj2;
bold1 = inpParse.Results.bold1;
bold2 = inpParse.Results.bold2;
runLabels = inpParse.Results.runLabels;
lambda = inpParse.Results.lambda;
lambdaVector = inpParse.Results.lambdaVector;
numRestarts = inpParse.Results.numRestarts;
zscore = inpParse.Results.zscore;
trustTheta = inpParse.Results.trustTheta;
trustDelta = inpParse.Results.trustDelta;
stepEpsilon = inpParse.Results.stepEpsilon;
quiet = inpParse.Results.quiet;
matchEpsilon = inpParse.Results.matchEpsilon;
allowNegative = inpParse.Results.allowNegative;

isUnset = @(y) any(cellfun(@(x) strcmp(x,y),inpParse.UsingDefaults));
if (isUnset('type'))
    error('type must be set to ''one'' or ''both'' ')
end
    
if (strcmp(type,'one'))
    bothRegionsFlag = 0;
    
    if (isUnset('adj1') || isUnset('bold1') || isUnset('bold2'))
        error('adj1, bold1, and bold2 required for type=one');
    elseif (~isUnset('numRestarts') || ~isUnset('trustTheta') || ...
            ~isUnset('trustDelta') || ~isUnset('stepEpsilon') || ...
            ~isUnset('adj2'))
        error(['numRestarts, trustTheta, trustDelta, adj2, and ' ...
            'stepEpsilon not applicable for type=one']);
    end
    
    if (isUnset('lambda'))
        if (isUnset('runLabels'))
            error(['either runLabels or lambda must be specified ' ...
                'for type=one']);
        end
        validateLambdaFlag = 1;
    else
        validateLambdaFlag = 0;
    end
    
elseif (strcmp(type,'both'))
    bothRegionsFlag = 1;
    
    if (isUnset('adj1') || isUnset('adj2') || isUnset('bold1') || ...
            isUnset('bold2') || ...
            isUnset('lambda'))
        error(['adj1, adj2, bold1, bold2, and lambda ' ...
            'required for type=both']);
    elseif (~isUnset('lambdaVector') || ~isUnset('runLabels'))
        error('lambdaVector and runLabels cannot be used with type=both');
    elseif (~isUnset('allowNegative'))
        error('Nonnegative constraint not implemented for ''both''');
    end
end
    
    
r1Size = size(adj1,1);
numTimepoints = size(bold1,2);

if (size(adj1,2) ~= r1Size)
    error('adj1 not square');
elseif (any(diag(adj1)))
    error('adj1 has diagonal elements');
elseif (any(any(~ismember(adj1,[0 1]))))
    error('adj1 has elements not equal to 0 or 1');
elseif (size(bold1,1) ~= r1Size)
    error('number of rows in bold1 must equal side length of adj1');
elseif (zscore ~= 0 && zscore ~= 1)
    error('zscore must be either 0 or 1');
elseif (quiet ~= 0 && quiet ~= 1)
    error('quiet must be either 0 or 1');
elseif (~isscalar(lambda) || lambda < 0)
    error('lambda must be a nonnegative scalar');
end

if (~bothRegionsFlag)
    if (validateLambdaFlag)
        if (~isrow(lambdaVector) || size(lambdaVector,2) < 1)
            error('lambdaVector must be a row vector');
        elseif (length(runLabels) ~= numTimepoints)
            error(['number of columns in bold1 must equal ' ...
                'length of runLabels']);
        elseif (~isrow(runLabels))
            error('runLabels must be a 1 x numTimepoints row vector');
        elseif (allowNegative ~= 0 && allowNegative ~= 1)
            error('allowNegative must be either 0 or 1');
        end
    end
else
    r2Size = size(adj2,1);
    if (size(adj2,2) ~= r2Size)
        error('adj2 not square');
    elseif (any(diag(adj2)))
        error('adj2 has diagonal elements');
    elseif (any(any(~ismember(adj2,[0 1]))))
        error('adj2 has elements not equal to 0 or 1');
    elseif (size(bold2,1) ~= r2Size)
        error('number of rows in bold2 must equal side length of adj2');
    elseif (~isscalar(numRestarts) || numRestarts <= 0)
        error('numRestarts must be a positive scalar');
    elseif (~isscalar(trustTheta) || trustTheta <= 0)
        error('trustTheta must be a positive scalar');
    elseif (~isscalar(trustDelta) || trustDelta <= 0)
        error('trustDelta must be a positive scalar');
    elseif (~isscalar(stepEpsilon) || stepEpsilon <= 0)
        error('stepEpsilon must be a positive scalar');
    elseif (~isscalar(matchEpsilon) || matchEpsilon <= 0)
        error('matchEpsilon must be a positive scalar');
    elseif (size(bold2,2) ~= numTimepoints)
        error('bold1 and bold2 must have the same number of columns');
    end
end



if (zscore)
    % Scale voxels by standard deviation
    
    for v1=1:size(bold1,1)
        bold1(v1,:) = (bold1(v1,:)-mean(bold1(v1,:)))/std(bold1(v1,:));
    end
    for v2=1:size(bold2,1)
        bold2(v2,:) = (bold2(v2,:)-mean(bold2(v2,:)))/std(bold2(v2,:));
    end
end


if (~bothRegionsFlag)
    % Learn connectivity map over one region
    
    spatialDiff = buildSpatialDiff(adj1,quiet);
    spatialDiff_aug = [spatialDiff zeros(size(spatialDiff,1),1)];
    bold1_aug = [bold1;ones(1,size(bold1,2))];
    
    if (validateLambdaFlag)
        % Tune lambda to maximize variance explained in cross-validation
        
        dispQuiet(['Finding best lambda out of ' ...
            num2str(length(lambdaVector)) ' possibilities:'],quiet);
        numRuns = length(unique(runLabels));
        bold2VarPerRun = norm(mean(mean(bold2,1)) - mean(bold2,1),2)^2/...
            numRuns;
        testFracExp = zeros(length(lambdaVector),1);
        for lambdaInd = 1:length(lambdaVector)
            valLambda = lambdaVector(lambdaInd);
            dispQuiet(['    Trying lambda = ' num2str(valLambda)],quiet);
            for leaveOutRun = 1:numRuns
                trainTRs = runLabels~=leaveOutRun;
                testTRs = runLabels==leaveOutRun;

                [a, b] = SingleRegionWeights(bold1_aug, bold2, ...
                                             spatialDiff_aug, ...
                                             valLambda, ...
                                             allowNegative, trainTRs);
                testFracExp(lambdaInd) = testFracExp(lambdaInd) + (1-...
                    norm(a'*bold1(:,testTRs)+b- ...
                    mean(bold2(:,testTRs),1),2)^2 ...
                    /(bold2VarPerRun))/numRuns;
            end
        end

        [~, maxLambdaInd] = max(testFracExp);
        lambda = lambdaVector(maxLambdaInd);
        dispQuiet(['Using lambda = ' num2str(lambda)],quiet);
    end
    
    % Estimate final connectivity map   
    [a, b] = SingleRegionWeights(bold1_aug, bold2, spatialDiff_aug, ...
                                 lambda, allowNegative, trainTRs);
    
    connWeights = a;
    fracVar = 1-norm(a'*bold1+b-mean(bold2,1),2)^2/...
        norm(mean(mean(bold2,1)) - mean(bold2,1),2)^2;
    dispQuiet(['Fraction of Variance explained: ' num2str(fracVar)],quiet);
    
else
    %Learn connectivity maps over both regions
    
    spatialDiff1 = buildSpatialDiff(adj1,quiet);
    spatialDiff2 = buildSpatialDiff(adj2,quiet);
    
    X = [bold1'./sqrt(numTimepoints) -1*bold2'./sqrt(numTimepoints);...
       sqrt(lambda)*spatialDiff1 zeros(size(spatialDiff1,1),r2Size);...
       zeros(size(spatialDiff2,1),r1Size) sqrt(lambda)*spatialDiff2];
   
   finalBetas = zeros(r1Size+r2Size,numRestarts);
   finalObj = zeros(1,numRestarts);
   for init=1:numRestarts
       % Try multiple initializations to find multiple solutions
       
       dispQuiet(['Starting initialization ' num2str(init)],quiet);
       beta = zeros(r1Size+r2Size,1);
       beta(randi([1 r1Size+r2Size])) = 1;
       
       while (1)
           lastbeta = beta;
           cvx_begin quiet
                cvx_solver sedumi
                variable s(r1Size+r2Size)
                minimize (norm(X*(lastbeta+s),2))
                subject to
                    norm(s,2) <= trustDelta;
                    abs((norm(lastbeta,2)^2-1) + 2*lastbeta'*s) ...
                        <= trustTheta;
                    (lastbeta+s) >= 0;
           cvx_end
           if(~strcmp(cvx_status,'Solved'))
               dispQuiet('Warning: optimization failure in cvx',quiet);
               break;
           end
           beta = lastbeta+s;
           
           dispQuiet(['    Obj = ' num2str(norm(X*(beta),2)^2)],quiet);
           if (norm(s,2) < stepEpsilon)
               break;
           end
       end
       
       finalBetas(:,init) = beta;
       finalObj(init) = norm(X*(beta),2)^2;
   end
   
   % Remove redundant solutions
   [~,sortedInds] = sort(finalObj,'ascend');
   uniqueInds = sortedInds([1 ...
       find(diff(finalObj(sortedInds))>matchEpsilon)+1]);
   
   connWeights = cell(2,1);
   connWeights{1} = finalBetas(1:r1Size,uniqueInds);
   connWeights{2} = finalBetas(r1Size+1:r1Size+r2Size,uniqueInds);
end


end

% Learn weights over region 1 to predict mean of region 2, using either
% real-valued weights or non-negative weights
function [a, b] = SingleRegionWeights(bold1_aug, bold2, spatialDiff_aug,...
                                      lambda, allowNegative, trainTRs)
    T = sum(trainTRs);
    target = mean(bold2(:,trainTRs),1)';
    C = [sparse(bold1_aug(:,trainTRs)'/sqrt(T)); ...
         sqrt(lambda)*spatialDiff_aug];
    H = C'*C;
    if (allowNegative)
        ab = H \ (bold1_aug(:,trainTRs)/T*target);
    else
        bound = sparse([zeros(size(bold1_aug,1)-1,1);-Inf]);
        opts_IPC = optimoptions('quadprog','Algorithm', ...
                   'interior-point-convex','display','off');
        weights_IPC = quadprog(H, ...
                sparse(-1*bold1_aug(:,trainTRs)/T*target), ...
                [], [], [], [], bound,[],[],opts_IPC);
        opts_TRR = optimoptions('quadprog','Algorithm', ...
                   'trust-region-reflective','display','off');
        ab = quadprog(H, ...
                sparse(-1*bold1_aug(:,trainTRs)/T*target), ...
                [], [], [], [], bound,[],weights_IPC,opts_TRR);
    end
    a = ab(1:(end-1));
    b = ab(end);
end

%Build a sparse connectivity matrix from adjacency matrix
function spDiff = buildSpatialDiff(adj,quiet)
    unconnectedFlag = 0;
    spDiff = zeros(sum(sum(adj)),size(adj,1));
    spInd = 1;
    for i = 1:size(adj,1)
        neighbors = find(adj(i,:));
        if (isempty(neighbors))
            unconnectedFlag = 1;
            continue;
        end
        nInd = sub2ind(size(spDiff), ...
            spInd:(spInd+length(neighbors)-1), neighbors);
        spDiff(nInd) = 1/sqrt(length(neighbors));
        spDiff(spInd:(spInd+length(neighbors)-1),i) = ...
            -1/sqrt(length(neighbors));
        spInd = spInd+length(neighbors);
    end
    spDiff = sparse(spDiff);
    if (unconnectedFlag)
        dispQuiet('Warning: some voxels have no neighbors',quiet);
    end
end


function dispQuiet(str,quietFlag)
    if (~quietFlag)
        disp(str);
    end
end