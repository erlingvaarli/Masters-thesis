% Cleaning data and plotting macrotexture

%% Import data
clear all

cs = {'notIndexed', crystalSymmetry('m-3m', [4.04 4.04 4.04], 'mineral', 'al')};
ssO = specimenSymmetry('orthorhombic');

datadir = '/Volumes/TOSHIBA EXT/mai 7/5000sek-12/5000sek-12-500X/1';
file = fullfile(datadir, 'Pattern.osc');

label = 'SUPATEST'; %HUSK Å ENDRE FOR NY PRØVE! 

ebsd = EBSD.load(file, cs, 'convertEuler2SpatialReferenceFrame','setting 2');
testebsd = ebsd;

% Set reference frame
setMTEXpref('xAxisDirection', 'north');
setMTEXpref('zAxisDirection', 'outOfPlane');
setMTEXpref('figSize','large')
% X || north        || ND
% Y || west         || RD
% Z || outOfPlane   || TD

%% Try and fail filtration
setMTEXpref('figSize','large')

iqfilterVal = 950; %Insert appropriate filterval
testebsd(testebsd.imagequality < iqfilterVal).phase = -1;
% 
figure
plot(testebsd)
%% Plotting the imagequality

figure
plot(ebsd,ebsd.imagequality)
%% If you need to reset filtration
testebsd = ebsd; 

%% Final filtration
ebsd(ebsd.imagequality < iqfilterVal).phase = -1;
%% Double check 
% PLotting OM 
setMTEXpref('FontSize',30)
setMTEXpref('figSize','huge')

oM2 = ipfHSVKey(ebsd('al'));
oM2.inversePoleFigureDirection = yvector;

br = orientation.byEuler(35*degree, 45*degree, 90*degree, cs{2}, ssO);
cu = orientation.byEuler(90*degree, 35*degree, 45*degree, cs, ssO);
cube = orientation.byEuler(0*degree, 0*degree, 0*degree, cs, ssO);
cubeND = orientation.byEuler(22*degree,0*degree,0*degree, cs, ssO);
cubeND45 = orientation.byEuler(45*degree, 0*degree, 0*degree, cs, ssO);
goss = orientation.byEuler(0, 45*degree, 0, cs, ssO);
p = orientation.byMiller([0 1 1], [1 2 2], cs, ssO);
q = orientation.byMiller([0 1 3], [2 3 1], cs, ssO);
s = orientation.byEuler(59*degree, 37*degree, 63*degree, cs, ssO);
componentss = {br,cu,cube,cubeND45,goss,p,s};


r = vector3d.Y;
components = [...
  orientation.brass(cs{2},ssO),...
  orientation.copper(cs{2},ssO),...
  orientation.cube(cs{2},ssO),...  
  orientation.byEuler(22*degree,0*degree,0*degree,cs{2},ssO),...
  orientation.goss(cs{2},ssO),...
  orientation.byEuler(70*degree,45*degree,0*degree,cs{2},ssO),...
  orientation.byEuler(59*degree,37*degree,63*degree,cs{2},ssO),...
  ];

markers = [...
    's',...
    '^',...
    'o',...
    'd',...
    'v',...
    'p',...
    '>',...
    ];

markercolors = [...
    'g',...
    'b',...
    'r',...
    'c',...
    'k',...
    'w',...
    'm'...
    ];
labels = {'Br','Cu','Cube','CubeND','Goss','P','S'};

figure
plot(ebsd('al'), oM2.orientation2color(ebsd('al').orientations))
%%
% ODF PART
cs = {'notIndexed', crystalSymmetry('m-3m', [4.04 4.04 4.04], 'mineral', 'al')};
ss = specimenSymmetry('1'); % Triclinic
ssO = specimenSymmetry('orthorhombic');

setMTEXpref('FontSize',18)
setMTEXpref('defaultColorMap','white2blackColorMap');
setMTEXpref('markerSize',12)
setMTEXpref('figSize','medium')

br = orientation.byEuler(35*degree, 45*degree, 90*degree, cs, ssO);
cu = orientation.byEuler(90*degree, 35*degree, 45*degree, cs, ssO);
cube = orientation.byEuler(0.001*degree, 0.001*degree, 0.001*degree, cs, ssO);
cubeND = orientation.byEuler(22*degree,0.001*degree,0.001*degree, cs, ssO);
cubeND45 = orientation.byEuler(45*degree, 0.001*degree, 0.001*degree, cs, ssO);
goss = orientation.byEuler(0, 45*degree, 0, cs, ssO);
p = orientation.byMiller([0 1 1], [1 2 2], cs, ssO);
q = orientation.byMiller([0 1 3], [2 3 1], cs, ssO);
s = orientation.byEuler(59*degree, 37*degree, 63*degree, cs, ssO);
idealoris = [br,cu,cube,goss,cubeND,p,s];

levelsODF = [0, 2, 4, 6, 8, 10, 15, 20, 25, 30];

% Calculate ODF
h = [Miller(1, 1, 1, cs{2}), Miller(1, 0, 0, cs{2}), Miller(1, 1, 0, cs{2})];
odf = calcDensity(ebsd('indexed').orientations, 'resolution', 7.5*degree,...
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
%
% Pole figures
pfAnnotations = @(varargin) text([vector3d.X, vector3d.Y],...
    {'RD', 'TD'}, 'BackgroundColor', 'w', 'tag', 'axesLabels', varargin{:});
setMTEXpref('pfAnnotations', pfAnnotations);
setMTEXpref('defaultColorMap','white2blackColorMap');
%% 
levelsPDF = [0,1,2,4,6,8];
figure
plotPDF(odf3, h, 'upper', 'projection', 'eangle', 'contourf',levelsPDF)
CLim(gcm,[0 8]);
print(fullfile('/Users/erlingaaresvarli/Library/Mobile Documents/com~apple~CloudDocs/Masteroppgave/OneDrive - NTNU/MATLAB/Resultater - MATLAB/Pole figures/Macro',append(label,'-polefig')),'-dpng');
% Exporting ODF and grain data
odfresultpath = '/Users/erlingaaresvarli/Library/Mobile Documents/com~apple~CloudDocs/Masteroppgave/OneDrive - NTNU/MATLAB/Resultater - MATLAB/ODF data/Macro';
fnameodf = fullfile(odfresultpath, append(label,'-ODFdata.mat'));
save(fnameodf,'odf')
% phi2 sections
setMTEXpref('FontSize',25)
setMTEXpref('figSize','medium')
setMTEXpref('defaultColorMap','white2blackColorMap');
setMTEXpref('markerSize',12)

levelsODF = [0,1,2,4,8,12,15,20,25,30];
limitODF = [0,10];

odf3.SS = specimenSymmetry('orthorhombic');
figure
plot(odf3, 'phi2', [0 45 65]*degree, 'contourf', levelsODF, 'minmax')
CLim(gcm,[0 30]);
%mtexColorbar
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

%% Fibres

% Writing fibres to file
f = fibre(cu, br, cs, ssO);
fibrepath_beta = '/Users/erlingaaresvarli/Library/Mobile Documents/com~apple~CloudDocs/Masteroppgave/OneDrive - NTNU/MATLAB/Resultater - MATLAB/Fibredata/Beta/Macro';

% generate list from fibres and evalute ODF at specific orientations
fibreOris = f.orientation;
evalOris = [];
evalIndex = [1 84 167 254 346 446 556 680 824 1000];
evalValues = zeros(1, 10);
for i=1:10
    ori = fibreOris(evalIndex(i));
    evalOris = [evalOris ori];
    evalValues(i) = eval(odf3, ori);
end

figure
plot(evalOris.phi2/degree, evalValues, '-o')
xlabel('\phi_2 \rightarrow', 'interpreter', 'tex')
ylabel('Orientation density f(g)', 'interpreter', 'tex')
xlim([45 90])

% Write fibre data to a csv file for further analysis
datafname = fullfile(fibrepath_beta, append(label,'-data_fibre_beta.csv'));

% Write header to file
fid = fopen(datafname, 'w');
fprintf(fid, '%s\r\n', 'phi1,Phi,phi2,fibreValue');
fclose(fid);

% Write Euler angles and intensities to file
dlmwrite(datafname, [(evalOris.phi1/degree)' (evalOris.Phi/degree)'...
    (evalOris.phi2/degree)' evalValues'], '-append')

%Plot intensity along fibre from Cube to Goss and write results to file
odf.SS = ssO;
f = fibre(cube, goss, cs, ssO);
fibrepath_cubegoss = '/Users/erlingaaresvarli/Library/Mobile Documents/com~apple~CloudDocs/Masteroppgave/OneDrive - NTNU/MATLAB/Resultater - MATLAB/Fibredata/Cube-Goss/Macro';


% Generate list from fibres and evalute ODF at specific orientations
fibreOris = f.orientation;
evalOris = [];
evalIndex = [1 111 222 333 444 555 666 777 888 1000];
evalValues = zeros(1, 10);
for i=1:10
    ori = fibreOris(evalIndex(i));
    evalOris = [evalOris ori];
    evalValues(i) = eval(odf3, ori);
end

figure
plot(evalOris.Phi/degree, evalValues, '-o')
xlabel('\Phi \rightarrow', 'interpreter', 'tex')
ylabel('Orientation density f(g)', 'interpreter', 'tex')
xlim([0 45])

% Write fibre data to a csv file for further analysis
datafname = fullfile(fibrepath_cubegoss, append(label,'data_fibre_cube_goss.csv'));


% Write header to file
fid = fopen(datafname, 'w');
fprintf(fid, '%s\r\n', 'phi1,Phi,phi2,fibreValue');
fclose(fid);

%Write Euler angles and intensities to file
dlmwrite(datafname, [(evalOris.phi1/degree)' (evalOris.Phi/degree)'...
    (evalOris.phi2/degree)' evalValues'], '-append')

%Plot intensity along fibre from Cube to ND-rotated Cube
odf.SS = ssO;
f = fibre(cube, cubeND45, cs, ssO);
fibrepath_cubecubeND45 = '/Users/erlingaaresvarli/Library/Mobile Documents/com~apple~CloudDocs/Masteroppgave/OneDrive - NTNU/MATLAB/Resultater - MATLAB/Fibredata/Cube-CubeND45/Macro';

 % Generate list from fibres and evalute ODF at specific orientations
fibreOris = f.orientation;
evalOris = [];
evalIndex = [1 111 222 333 444 555 666 777 888 1000];
evalValues = zeros(1,10);
for i=1:10
    ori = fibreOris(evalIndex(i));
    evalOris = [evalOris ori];
    evalValues(i) = eval(odf3, ori);
end

figure
plot(evalOris.phi1/degree, evalValues, '-o')
xlabel('\phi_1 \rightarrow', 'interpreter', 'tex')
ylabel('Orientation density f(g)', 'interpreter', 'tex')
xlim([0 45])

%Write fibre data to csv file for further analysis in Python
datafname = fullfile(fibrepath_cubecubeND45, append(label,'data_fibre_cube_cubeND45.csv'));

fid = fopen(datafname, 'w');
fprintf(fid, '%s\r\n', 'phi1,Phi,phi2,fibreValue');
fclose(fid);

dlmwrite(datafname, [(evalOris.phi1/degree)' (evalOris.Phi/degree)'...
    (evalOris.phi2/degree)' evalValues'], '-append')
