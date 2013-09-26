function connectivityDemo()
%connectivityDemo Gives two examples of using learnConnectivity

if (~exist('cvx_version'))
    error('learnConnectivity requires CVX - see http://cvxr.com/cvx/');
end


% Example One: Learn a connectivity map over hV4, using the mean PPA
% timecourse as the seed. The result shows that hV4 voxels which respond to
% higher eccentricities are more strongly connected to PPA than those
% voxels whose receptive fields are more foveal. The lambda value used in
% the optimization is automatically selected based on the fraction of PPA
% variance explained during cross-validation.

load('sampleDataV4PPA.mat');

% Construct adjacency matrix, using nearest neighbors on cortical flatmap
N = 10;
numV4voxels = size(V4bold,1);
adjV4 = zeros(numV4voxels,numV4voxels);
for vox = 1:numV4voxels
    [~, closestInds] = sort(norms( ...
        V4FlatmapLocations(vox*ones(numV4voxels,1),:)-V4FlatmapLocations...
        ,2,2), 'ascend');
    
    adjV4(vox, closestInds(2:N+1)) = 1;
end


% Learn a connectivity map over hV4, for connectivity to PPA
V4Weights = learnConnectivity('type','one','adj1',adjV4,'bold1',V4bold,...
    'bold2',PPAbold,'runLabels',runLabels,'lambdaVector',logspace(-4,5,10));

% Draws a figure showing hV4 eccentricity and the learned weights
figure('Position',[200 800 1200 500]);
subplot(1,2,1); hold on; axis([-150 125 -125 -30]);
drawConvHull([V1FlatmapLocations(V1FlatmapLocations(:,1)<0,:); ...
    V2FlatmapLocations(V2FlatmapLocations(:,1)<0,:); ...
    VPFlatmapLocations(VPFlatmapLocations(:,1)<0,:)]);
text(-40,-70,'lV1,lV2,lVP');
drawConvHull([V1FlatmapLocations(V1FlatmapLocations(:,1)>0,:); ...
    V2FlatmapLocations(V2FlatmapLocations(:,1)>0,:); ...
    VPFlatmapLocations(VPFlatmapLocations(:,1)>0,:)]);
text(-15,-90,'rV1,rV2,rVP');
drawConvHull(V4FlatmapLocations(V4FlatmapLocations(:,1)<0,:));
text(-130,-90,'lV4');
drawConvHull(V4FlatmapLocations(V4FlatmapLocations(:,1)>0,:));
text(95,-90,'rV4');
scatter(V4FlatmapLocations(:,1),V4FlatmapLocations(:,2),40,...
    V4Eccentricity,'filled');
title('Eccentricity of V4 Receptive Fields');
colorbar;

subplot(1,2,2); hold on; axis([-150 125 -125 -30]);
[~,weightInds] = sort(V4Weights,'descend');
drawConvHull([V1FlatmapLocations(V1FlatmapLocations(:,1)<0,:); ...
    V2FlatmapLocations(V2FlatmapLocations(:,1)<0,:); ...
    VPFlatmapLocations(VPFlatmapLocations(:,1)<0,:)]);
text(-40,-70,'lV1,lV2,lVP');
drawConvHull([V1FlatmapLocations(V1FlatmapLocations(:,1)>0,:); ...
    V2FlatmapLocations(V2FlatmapLocations(:,1)>0,:); ...
    VPFlatmapLocations(VPFlatmapLocations(:,1)>0,:)]);
text(-15,-90,'rV1,rV2,rVP');
drawConvHull(V4FlatmapLocations(V4FlatmapLocations(:,1)<0,:));
text(-130,-90,'lV4');
drawConvHull(V4FlatmapLocations(V4FlatmapLocations(:,1)>0,:));
text(95,-90,'rV4');
scatter(V4FlatmapLocations(weightInds,1),V4FlatmapLocations(weightInds,2),...
    40,V4Weights(weightInds),'filled');
title('Connectivity Weights with PPA');
colorbar;





% Example Two: Learn connectivity maps over V1 and VP (ventral V3)
% simultaneously, to discover which sets of voxels in each region are most
% connected to each other. The random initialization is seeded in this
% example so that the first two initializations return two different
% solutions, one in the left hemisphere and one in the right hemisphere.

load('sampleDataV1VP.mat');

% Construct adjacency matrices, using nearest neighbors on cortical flatmap
N = 10;
numV1voxels = size(V1bold,1);
adjV1 = zeros(numV1voxels,numV1voxels);
for vox = 1:numV1voxels
    [~, closestInds] = sort(norms( ...
        V1FlatmapLocations(vox*ones(numV1voxels,1),:)-V1FlatmapLocations...
        ,2,2), 'ascend');
    
    adjV1(vox, closestInds(2:N+1)) = 1;
end
N = 10;
numVPvoxels = size(VPbold,1);
adjVP = zeros(numVPvoxels,numVPvoxels);
for vox = 1:numVPvoxels
    [~, closestInds] = sort(norms( ...
        VPFlatmapLocations(vox*ones(numVPvoxels,1),:)-VPFlatmapLocations...
        ,2,2), 'ascend');
    
    adjVP(vox, closestInds(2:N+1)) = 1;
end

% Set random seed such that only two initializations (restarts) are needed
% to find two different solutions
matlabVer = version('-release');
if (str2num(matlabVer(1:4)) <= 2010)
    RandStream.setDefaultStream(RandStream('mt19937ar','Seed',2));
else
    RandStream.setGlobalStream(RandStream('mt19937ar','Seed',2));
end

% Learn connectivity maps over V1 and VP
connWeights = learnConnectivity('type','both','adj1',adjV1,'bold1',V1bold,...
    'adj2',adjVP,'bold2',VPbold,'lambda',100,'zscore',1,'numRestarts',2);

% Draw a figure showing both solutions
figure('Position',[400 600 1200 500]);
subplot(1,2,1); hold on;
scatter(V1FlatmapLocations(:,1),V1FlatmapLocations(:,2),30,...
    connWeights{1}(:,1),'filled');
scatter(VPFlatmapLocations(:,1),VPFlatmapLocations(:,2),30,...
    connWeights{2}(:,1),'filled');
axis([-150 125 0 120]); colorbar;
title('Solution 1');
drawConvHull(V1FlatmapLocations(V1FlatmapLocations(:,1)<0,:));
text(-40,60,'lV1');
drawConvHull(V1FlatmapLocations(V1FlatmapLocations(:,1)>0,:));
text(10,60,'rV1');
drawConvHull(VPFlatmapLocations(VPFlatmapLocations(:,1)<0,:));
text(-110,30,'lVP');
drawConvHull(VPFlatmapLocations(VPFlatmapLocations(:,1)>0,:));
text(100,30,'lVP');


subplot(1,2,2); hold on;
scatter(V1FlatmapLocations(:,1),V1FlatmapLocations(:,2),30,...
    connWeights{1}(:,2),'filled');
scatter(VPFlatmapLocations(:,1),VPFlatmapLocations(:,2),30,...
    connWeights{2}(:,2),'filled');
axis([-150 125 0 120]); colorbar;
title('Solution 2');
drawConvHull(V1FlatmapLocations(V1FlatmapLocations(:,1)<0,:));
text(-40,60,'lV1');
drawConvHull(V1FlatmapLocations(V1FlatmapLocations(:,1)>0,:));
text(10,60,'rV1');
drawConvHull(VPFlatmapLocations(VPFlatmapLocations(:,1)<0,:));
text(-110,30,'lVP');
drawConvHull(VPFlatmapLocations(VPFlatmapLocations(:,1)>0,:));
text(100,30,'lVP');
end

% Utility function for outlining a region of interest
function drawConvHull(pts)
    outline = convhull(pts(:,1),pts(:,2));
    h = line(pts(outline,1),pts(outline,2));
    set(h,'Color','k');
    set(h,'LineWidth',2);
end
