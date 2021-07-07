%% Import data
clear all

% Defining crystal system (cs) and specimen symmetry (ssO)
cs = {'notIndexed', crystalSymmetry('m-3m', [4.04 4.04 4.04], 'mineral', 'al')};
ssO = specimenSymmetry('orthorhombic');

% Input file
datadir = '/Users/erlingaaresvarli/Desktop/HD Micro OSC files';
file = fullfile(datadir, '3 - 784713-Midt-16mars-Pattern_sda.osc');

% Loading EBSD data
ebsd = EBSD.load(file, cs, 'convertEuler2SpatialReferenceFrame','setting 2');
testebsd = ebsd;

% Set reference frame
setMTEXpref('xAxisDirection', 'north');
setMTEXpref('zAxisDirection', 'outOfPlane');
% X || north        || ND
% Y || west         || RD
% Z || outOfPlane   || TD


%% Filter away low IQ-measurements (Second phase particle data)

iqfilterVal = 2000; % Insert appropriate filterval, compare result with SE image so that resulting notIndex areas correspond with second phase particles
ebsd(ebsd.imagequality < iqfilterVal).phase = -1;
%% Grain reconstruction
% Reconstructing grains, using a treshold of 2 degrees
[grains,ebsd.grainId] = calcGrains(ebsd,'angle',2*degree);

if isempty(grains('notIndexed')) == 0
    %Removing notIndexed grains with a boundary pixels/total pixel ratio > 0.25
    notIndexed = grains('notIndexed');
    toRemove = notIndexed((log(notIndexed.grainSize ./ notIndexed.boundarySize))<-0.6);
    ebsd(toRemove) = [];
    
    % Reconstructing grains, using a treshold of 2 degrees 
    [grains,ebsd.grainId,ebsd.mis2mean] = calcGrains(ebsd,'angle',2*degree);
end

%Remove grains with less than 10 pixels
valid_grains = grains(grains.grainSize >= 10);
ebsd_clean = ebsd(valid_grains);
[newgrains,ebsd_clean.grainId,ebsd_clean.mis2mean] = calcGrains(ebsd_clean,'angle',2*degree);
%% Orientation map 
% Defining inverse polefigure key
oM2 = ipfHSVKey(ebsd_clean('al'));
oM2.inversePoleFigureDirection = yvector;

% Plotting orientation map
figure
plot(ebsd_clean('al'), oM2.orientation2color(ebsd_clean('al').orientations))
hold on
plot(newgrains.boundary)
hold off

%% Grain segmentation
% Segmenting grains and plotting segmentation maps
[grains, grainsOR, grainsSub, grainsRex, graindata] = ebsd_fraction_recrystallized_with_og(newgrains);

% Calculate fraction recrystallized
Xrex = sum(grains(grains.RX == 1).area) / sum(grains.area);

% Calculating fraction of Subgrains
Xsub = sum(grains(grains.RX == 0).area) / sum(grains.area);

% Calculate area weighed grain sizes 
average_grain_size_Sub = 'NaN';
average_grain_size_RX = 'NaN';
average_grain_size_OG = 'NaN';

if isempty(grainsSub) == 0
    average_grain_size_Sub = sum(grainsSub.area .* grainsSub.prop.ECD)./sum(grainsSub.area);
end

if isempty(grainsRex) == 0
    average_grain_size_RX = sum(grainsRex.area .* grainsRex.prop.ECD)./sum(grainsRex.area);
end

if isempty(grainsOR) == 0
    average_grain_size_OG = sum(grainsOR.area .* grainsOR.prop.ECD)./sum(grainsOR.area);
end

%% Calculate ODF
% Setting crystal symmetry and speciment symmetry

cs = {'notIndexed', crystalSymmetry('m-3m', [4.04 4.04 4.04], 'mineral', 'al')};
ssO = specimenSymmetry('orthorhombic');

% Calculate ODF
h = [Miller(1, 1, 1, cs{2}), Miller(1, 0, 0, cs{2}), Miller(1, 1, 0, cs{2})];
odf = calcDensity(ebsd_clean('indexed').orientations, 'resolution', 7.5*degree,...
    'halfwidth', 5*degree);

% Rotate ODF so that these axes are parallel
% X || north        || RD: Old Y to new X: [0 1 0]
% Y || west         || TD: Old Z to new Y: [0 0 1]
% Z || OutOfPlane   || ND: Old X to new Z: [1 0 0]

rotODF = rotation.byMatrix([0 1 0; 0 0 1; 1 0 0]);
odf2 = rotate(odf, rotODF);

% Final rotation correcting for a slight misalignment of the sample in the
% microscope
[odf3, rot_inv3] = centerSpecimen(odf2, yvector);

%% Plotting phi2 sections
% Setting contour levels for the phi2 sections
levelsODF = [0, 2, 4, 6, 8, 10, 15, 20, 25, 30];
phi2_colorlim = [0 30];

% Setting correct specimen symmetry for the ODF
odf3.SS = specimenSymmetry('orthorhombic');

% Defining texture components
br = orientation.byEuler(35*degree, 45*degree, 90*degree, cs{2}, ssO);
cu = orientation.byEuler(90*degree, 35*degree, 45*degree, cs, ssO);
cube = orientation.byEuler(0*degree, 0*degree, 0*degree, cs, ssO);
cubeND = orientation.byEuler(22*degree,0*degree,0*degree, cs, ssO);
cubeND45 = orientation.byEuler(45*degree, 0*degree, 0*degree, cs, ssO);
goss = orientation.byEuler(0, 45*degree, 0, cs, ssO);
p = orientation.byMiller([0 1 1], [1 2 2], cs, ssO);
q = orientation.byMiller([0 1 3], [2 3 1], cs, ssO);
s = orientation.byEuler(59*degree, 37*degree, 63*degree, cs, ssO);

% Plotting phi2 sections and annotating texture components
figure
plot(odf3, 'phi2', [0 45 65]*degree, 'contourf', levelsODF, 'minmax')
CLim(gcm,phi2_colorlim);

hold on
annotate(br.symmetrise, 'marker', 'd', 'markerfacecolor', 'g','DisplayName','Br')
annotate(cu.symmetrise, 'marker', '^', 'markerfacecolor', 'b','DisplayName','Cu')
annotate(cube.symmetrise, 'marker', 's', 'markerfacecolor', 'r','DisplayName','Cube')
annotate(cubeND.symmetrise, 'marker', 's', 'markerfacecolor', 'c','DisplayName','CubeND')
annotate(goss.symmetrise, 'marker', 'o', 'markerfacecolor', 'k','DisplayName','Goss')
annotate(p.symmetrise, 'marker', '>', 'markerfacecolor',[1 0.6 0] ,'DisplayName','P')
annotate(s.symmetrise, 'marker', 'd', 'markerfacecolor', 'm','DisplayName','S')
hold on
[h,icons] = legend('Location','southoutside','Orientation','Horizontal');
n = ceil(numel(icons)/2);
newicons= icons(n+1:end);
for k=1:(length(newicons))
    newicons(k).Children.MarkerSize = 12;  
end
hold off
