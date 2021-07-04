% %Plot Microtexture data
% %Importing AVERAGE DOT PRODUCT MAP
% adp = readtable('/Volumes/TOSHIBA EXT/DOTPRODUCTMAP_TEST-2.csv','ReadVariableNames',1);
% 
% xcords = [0;0.0909999981522560;0.181999996304512;0.273000001907349;0.363999992609024;0.455000013113022;0.546000003814697;0.637000024318695;0.727999985218048;0.819000005722046;0.910000026226044;1.00100004673004;1.09200000762939;1.18299996852875;1.27400004863739;1.36500000953674;1.45599997043610;1.54700005054474;1.63800001144409;1.72899997234344;1.82000005245209;1.91100001335144;2.00200009346008;2.09299993515015;2.18400001525879;2.27500009536743;2.36599993705750;2.45700001716614;2.54800009727478;2.63899993896484;2.73000001907349;2.82100009918213;2.91199994087219;3.00300002098084;3.09400010108948;3.18499994277954;3.27600002288818;3.36700010299683;3.45799994468689;3.54900002479553;3.64000010490418;3.73099994659424;3.82200002670288;3.91300010681152;4.00400018692017;4.09499979019165;4.18599987030029;4.27699995040894;4.36800003051758;4.45900011062622;4.55000019073486;4.64099979400635;4.73199987411499;4.82299995422363;4.91400003433228;5.00500011444092;5.09600019454956;5.18699979782105;5.27799987792969;5.36899995803833;5.46000003814697;5.55100011825562;5.64200019836426;5.73299980163574;5.82399988174439;5.91499996185303;6.00600004196167;6.09700012207031;6.18800020217896;6.27899980545044;6.36999988555908;6.46099996566773;6.55200004577637;6.64300012588501;6.73400020599365;6.82499980926514;6.91599988937378;7.00699996948242;7.09800004959106;7.18900012969971;7.28000020980835;7.37099981307983;7.46199989318848;7.55299997329712;7.64400005340576;7.73500013351440;7.82600021362305;7.91699981689453;8.00800037384033;8.09899997711182;8.18999958038330;8.28100013732910;8.37199974060059;8.46300029754639;8.55399990081787;8.64500045776367;8.73600006103516;8.82699966430664;8.91800022125244;9.00899982452393;9.10000038146973;9.19099998474121;9.28199958801270;9.37300014495850;9.46399974822998;9.55500030517578;9.64599990844727;9.73700046539307;9.82800006866455;9.91899967193604;10.0100002288818;10.1009998321533;10.1920003890991;10.2829999923706;10.3739995956421;10.4650001525879;10.5559997558594;10.6470003128052;10.7379999160767;10.8290004730225;10.9200000762939;11.0109996795654;11.1020002365112;11.1929998397827;11.2840003967285;11.3750000000000;11.4659996032715;11.5570001602173;11.6479997634888;11.7390003204346;11.8299999237061;11.9209995269775;12.0120000839233;12.1029996871948;12.1940002441406;12.2849998474121;12.3760004043579;12.4670000076294;12.5579996109009;12.6490001678467;12.7399997711182;12.8310003280640;12.9219999313355;13.0129995346069;13.1040000915527;13.1949996948242;13.2860002517700;13.3769998550415;13.4680004119873;13.5590000152588;13.6499996185303;13.7410001754761;13.8319997787476;13.9230003356934;14.0139999389648;14.1049995422363;14.1960000991821;14.2869997024536;14.3780002593994;14.4689998626709;14.5600004196167;14.6510000228882;14.7419996261597;14.8330001831055;14.9239997863770;15.0150003433228;15.1059999465942;15.1969995498657;15.2880001068115;15.3789997100830;15.4700002670288;15.5609998703003;15.6520004272461;15.7430000305176;15.8339996337891;15.9250001907349;16.0160007476807;16.1070003509522;16.1979999542236;16.2889995574951;16.3799991607666;16.4710006713867;16.5620002746582;16.6529998779297;16.7439994812012;16.8349990844727;16.9260005950928;17.0170001983643;17.1079998016357;17.1989994049072;17.2900009155273;17.3810005187988;17.4720001220703;17.5629997253418;17.6539993286133;17.7450008392334;17.8360004425049;17.9270000457764;18.0179996490479;18.1089992523193;18.2000007629395;18.2910003662109;18.3819999694824;18.4729995727539;18.5639991760254;18.6550006866455;18.7460002899170;18.8369998931885;18.9279994964600;19.0189990997314;19.1100006103516;19.2010002136230;19.2919998168945;19.3829994201660;19.4740009307861;19.5650005340576;19.6560001373291;19.7469997406006;19.8379993438721;19.9290008544922;20.0200004577637;20.1110000610352;20.2019996643066;20.2929992675781;20.3840007781982;20.4750003814697;20.5659999847412;20.6569995880127;20.7479991912842;20.8390007019043;20.9300003051758;21.0209999084473;21.1119995117188;21.2029991149902;21.2940006256104;21.3850002288818;21.4759998321533;21.5669994354248;21.6580009460449;21.7490005493164;21.8400001525879;21.9309997558594;22.0219993591309;22.1130008697510;22.2040004730225;22.2950000762939;22.3859996795654;22.4769992828369;22.5680007934570;22.6590003967285;22.7500000000000;22.8409996032715;22.9319992065430;23.0230007171631;23.1140003204346;23.2049999237061;23.2959995269775;23.3869991302490;23.4780006408691;23.5690002441406;23.6599998474121;23.7509994506836;23.8419990539551;23.9330005645752;24.0240001678467;24.1149997711182;24.2059993743897;24.2970008850098;24.3880004882813;24.4790000915527;24.5699996948242;24.6609992980957;24.7520008087158;24.8430004119873;24.9340000152588;25.0249996185303;25.1159992218018;25.2070007324219;25.2980003356934;25.3889999389648;25.4799995422363;25.5709991455078;25.6620006561279;25.7530002593994;25.8439998626709;25.9349994659424;26.0259990692139;26.1170005798340;26.2080001831055;26.2989997863770;26.3899993896484;26.4810009002686;26.5720005035400;26.6630001068115;26.7539997100830;26.8449993133545;26.9360008239746;27.0270004272461;27.1180000305176;27.2089996337891;27.2999992370605;27.3910007476807;27.4820003509522;27.5729999542236;27.6639995574951;27.7549991607666;27.8460006713867;27.9370002746582;28.0279998779297;28.1189994812012;28.2099990844727;28.3010005950928;28.3920001983643;28.4829998016357;28.5739994049072;28.6650009155273;28.7560005187988;28.8470001220703;28.9379997253418;29.0289993286133;29.1200008392334;29.2110004425049;29.3020000457764;29.3929996490479;29.4839992523193;29.5750007629395;29.6660003662109;29.7569999694824;29.8479995727539;29.9389991760254;30.0300006866455;30.1210002899170;30.2119998931885;30.3029994964600;30.3939990997314;30.4850006103516;30.5760002136230;30.6669998168945;30.7579994201660;30.8490009307861;30.9400005340576;31.0310001373291;31.1219997406006;31.2129993438721;31.3040008544922;31.3950004577637;31.4860000610352;31.5769996643066;31.6679992675781;31.7590007781982;31.8500003814697;31.9409999847412;32.0320014953613;32.1230010986328;32.2140007019043;32.3050003051758;32.3959999084473;32.4869995117188;32.5779991149902;32.6689987182617;32.7599983215332;32.8510017395020;32.9420013427734;33.0330009460449;33.1240005493164;33.2150001525879;33.3059997558594;33.3969993591309;33.4879989624023;33.5789985656738;33.6699981689453;33.7610015869141;33.8520011901856;33.9430007934570;34.0340003967285;34.1250000000000;34.2159996032715;34.3069992065430;34.3979988098145;34.4889984130859;34.5800018310547;34.6710014343262;34.7620010375977;34.8530006408691;34.9440002441406;35.0349998474121;35.1259994506836;35.2169990539551;35.3079986572266;35.3989982604981;35.4900016784668;35.5810012817383;35.6720008850098;35.7630004882813;35.8540000915527;35.9449996948242;36.0359992980957;36.1269989013672;36.2179985046387;36.3089981079102;36.4000015258789;36.4910011291504;36.5820007324219;36.6730003356934;36.7639999389648;36.8549995422363;36.9459991455078;37.0369987487793;37.1279983520508;37.2190017700195;37.3100013732910;37.4010009765625;37.4920005798340;37.5830001831055;37.6739997863770;37.7649993896484;37.8559989929199;37.9469985961914;38.0379981994629;38.1290016174316;38.2200012207031;38.3110008239746;38.4020004272461;38.4930000305176;38.5839996337891;38.6749992370606;38.7659988403320;38.8569984436035;38.9480018615723;39.0390014648438;39.1300010681152;39.2210006713867;39.3120002746582;39.4029998779297;39.4939994812012;39.5849990844727;39.6759986877441;39.7669982910156;39.8580017089844;39.9490013122559;40.0400009155273;40.1310005187988;40.2220001220703;40.3129997253418;40.4039993286133;40.4949989318848;40.5859985351563;40.6769981384277;40.7680015563965;40.8590011596680;40.9500007629395;41.0410003662109;41.1319999694824;41.2229995727539;41.3139991760254;41.4049987792969;41.4959983825684;41.5870018005371;41.6780014038086;41.7690010070801;41.8600006103516;41.9510002136231;42.0419998168945;42.1329994201660;42.2239990234375;42.3149986267090;42.4059982299805;42.4970016479492;42.5880012512207;42.6790008544922;42.7700004577637;42.8610000610352;42.9519996643066;43.0429992675781;43.1339988708496;43.2249984741211;43.3160018920898;43.4070014953613;43.4980010986328;43.5890007019043;43.6800003051758;43.7709999084473;43.8619995117188;43.9529991149902;44.0439987182617;44.1349983215332;44.2260017395020;44.3170013427734;44.4080009460449;44.4990005493164;44.5900001525879;44.6809997558594;44.7719993591309;44.8629989624023;44.9539985656738;45.0449981689453;45.1360015869141;45.2270011901856;45.3180007934570;45.4090003967285;45.5000000000000;45.5909996032715;45.6819992065430;45.7729988098145;45.8639984130859;45.9550018310547;46.0460014343262;46.1370010375977;46.2280006408691;46.3190002441406;46.4099998474121;46.5009994506836;46.5919990539551;46.6829986572266;46.7739982604981;46.8650016784668;46.9560012817383;47.0470008850098;47.1380004882813;47.2290000915527;47.3199996948242;47.4109992980957;47.5019989013672;47.5929985046387;47.6839981079102;47.7750015258789;47.8660011291504;47.9570007324219;48.0480003356934;48.1389999389648;48.2299995422363;48.3209991455078;48.4119987487793;48.5029983520508;48.5940017700195;48.6850013732910;48.7760009765625;48.8670005798340;48.9580001831055;49.0489997863770;49.1399993896484;49.2309989929199;49.3219985961914;49.4129981994629;49.5040016174316;49.5950012207031;49.6860008239746;49.7770004272461;49.8680000305176;49.9589996337891;50.0499992370606;50.1409988403320;50.2319984436035;50.3230018615723;50.4140014648438;50.5050010681152;50.5960006713867;50.6870002746582;50.7779998779297;50.8689994812012;50.9599990844727;51.0509986877441;51.1419982910156;51.2330017089844;51.3240013122559;51.4150009155273;51.5060005187988;51.5970001220703;51.6879997253418;51.7789993286133;51.8699989318848;51.9609985351563;52.0519981384277;52.1430015563965;52.2340011596680;52.3250007629395;52.4160003662109;52.5069999694824;52.5979995727539;52.6889991760254;52.7799987792969;52.8709983825684;52.9620018005371;53.0530014038086;53.1440010070801;53.2350006103516;53.3260002136231;53.4169998168945;53.5079994201660;53.5989990234375;53.6899986267090;53.7809982299805;53.8720016479492;53.9630012512207;54.0540008544922;54.1450004577637;54.2360000610352;54.3269996643066;54.4179992675781;54.5089988708496;54.5999984741211;54.6910018920898;54.7820014953613;54.8730010986328;54.9640007019043;55.0550003051758;55.1459999084473;55.2369995117188;55.3279991149902;55.4189987182617;55.5099983215332;55.6010017395020;55.6920013427734;55.7830009460449;55.8740005493164;55.9650001525879;56.0559997558594;56.1469993591309;56.2379989624023;56.3289985656738;56.4199981689453;56.5110015869141;56.6020011901856;56.6930007934570;56.7840003967285;56.8750000000000;56.9659996032715;57.0569992065430;57.1479988098145;57.2389984130859;57.3300018310547;57.4210014343262;57.5120010375977;57.6030006408691;57.6940002441406;57.7849998474121;57.8759994506836;57.9669990539551;58.0579986572266;58.1489982604981;58.2400016784668;58.3310012817383;58.4220008850098;58.5130004882813;58.6040000915527;58.6949996948242;58.7859992980957;58.8769989013672;58.9679985046387;59.0589981079102;59.1500015258789;59.2410011291504;59.3320007324219;59.4230003356934;59.5139999389648;59.6049995422363;59.6959991455078;59.7869987487793;59.8779983520508;59.9690017700195;60.0600013732910;60.1510009765625;60.2420005798340;60.3330001831055;60.4239997863770;60.5149993896484]
% 
% %%
% a = zeros(666,1);
% aa = 666*666
% %%
% b = zeros(1,443556);
% c = linspace(0.0, 60.5149993896484, 666);
% a(1) = 1;
% for i=2:666
%     a(i) = a(i-1)+666;
% end 
% 
% for i=1:666
%     b(1,a(i):a(i+1)-1) = c(:);
% end
% 

%% Import data
clear all

cs = {'notIndexed', crystalSymmetry('m-3m', [4.04 4.04 4.04], 'mineral', 'al')};
ssO = specimenSymmetry('orthorhombic');

datadir = '/Volumes/TOSHIBA EXT/april 23/4703-RT/2000X/1';
file = fullfile(datadir, 'april 23-4703-RT-2000X-Pattern_sda.osc');
label = 'april 23-4703-RT-2000X-ny'; %HUSK Å ENDRE FOR NY PRØVE! 


%REMEMBER%: Use length = 30 um in scalebar for RX-samples!
%REMEMBER%: Use length = 30 um in scalebar for RX-samples!
%REMEMBER%: Use length = 30 um in scalebar for RX-samples!
%REMEMBER%: Use length = 30 um in scalebar for RX-samples!
%REMEMBER%: Use length = 30 um in scalebar for RX-samples!


ebsd = EBSD.load(file, cs, 'convertEuler2SpatialReferenceFrame','setting 2');
testebsd = ebsd;
% Set reference frame
setMTEXpref('xAxisDirection', 'north');
setMTEXpref('zAxisDirection', 'outOfPlane');
setMTEXpref('figSize','large')
% X || north        || ND
% Y || west         || RD
% Z || outOfPlane   || TD

%% If wrong scalefactor is used:
ebsd.prop.x = ebsd.prop.x/1.65;
ebsd.prop.y = ebsd.prop.y/1.65;
ebsd.unitCell = ebsd.unitCell/1.65;

%%
% Filter away low IQ-measurements
iqfilterVal = 2000; %Insert appropriate filterval

testebsd(testebsd.imagequality  < iqfilterVal).phase = -1;


figure
plot(testebsd)
% testebsd(testebsd.fit < fitval).phase = -1;
% figure
% plot(testebsd)




%% If you need to reset filtration


testebsd = ebsd;
%testebsd.confidenceindex = testebsd.imagequality .* testebsd.confidenceindex;


%% Final filtration: 
ebsd(ebsd.imagequality < iqfilterVal).phase = -1;


% Calculating grains
[grains,ebsd.grainId] = calcGrains(ebsd,'angle',2*degree);

% Removing notIndexed grains with a certain 
if isempty(grains('notIndexed')) == 0
    %[grains,ebsd.grainId,ebsd.mis2mean] = calcGrains(ebsd); OLD
    notIndexed = grains('notIndexed');
    %
    toRemove = notIndexed((log(notIndexed.grainSize ./ notIndexed.boundarySize))<-0.6);

    %Se over denne!!
    ebsd(toRemove) = [];
    [grains,ebsd.grainId,ebsd.mis2mean] = calcGrains(ebsd,'angle',2*degree);
end
ebsd2 = ebsd;

 
%Remove grains with less than 10 pixels
valid_grains = grains(grains.grainSize >= 10);
ebsd_clean = ebsd(valid_grains);
[newgrains,ebsd_clean.grainId,ebsd_clean.mis2mean] = calcGrains(ebsd_clean,'angle',2*degree);

% PLotting OM and IPF with annotated texture components
%setMTEXpref('figSize','large')
setMTEXpref('FontSize',30)
setMTEXpref('figSize','huge')

oM2 = ipfHSVKey(ebsd_clean('al'));
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

% For plotting IPF KEY
% figure
% plot(oM2)
% hold on
% for i = 1:length(components)
%   plotIPDF(components(i),r,'Marker',markers(i),'MarkerSize',14,'MarkerEdgeColor','k','MarkerColor', markercolors(i),...
%     'DisplayName',labels{i})
%   hold on
% end
% hold on
% [h,icons] = legend('Location','eastoutside','Orientation','Vertical');
% n = ceil(numel(icons)/2);
% newicons= icons(n+1:end);
% for k=1:(length(newicons))
%     newicons(k).Children.MarkerSize = 14;  
% end
% hold off
% 
%testgrains = grains(newgrains.phase == 0);
%%
figure
plot(ebsd_clean('al'), oM2.orientation2color(ebsd_clean('al').orientations))
hold on
%plot(testgrains.boundary)
plot(newgrains.boundary)
hold off
%
%ompath = '/Users/erlingaaresvarli/Library/Mobile Documents/com~apple~CloudDocs/Masteroppgave/MATLAB/Resultater - MATLAB/Orientation maps/Back-Annealed/';
ompath = '/Users/erlingaaresvarli/Desktop/'
filenameom = append(ompath,label,'-om');
print(filenameom,'-dpdf','-bestfit')
 
% Finding ideal orientations

cs = {'notIndexed', crystalSymmetry('m-3m', [4.04 4.04 4.04], 'mineral', 'al')};
ssO = specimenSymmetry('orthorhombic');
grains = newgrains;


%
% Set reference frame
setMTEXpref('xAxisDirection', 'north');
setMTEXpref('zAxisDirection', 'outOfPlane');
setMTEXpref('figSize','large')
setMTEXpref('FontSize',18)

% X || north        || ND
% Y || west         || RD
% Z || outOfPlane   || TD


labels = {'Br','Cu','Cube','CubeND','Goss','P','S'};
idealOris = {...
  orientation.brass(cs{2},ssO),...
  orientation.copper(cs{2},ssO),...
  orientation.cube(cs{2},ssO),...  
  orientation.byEuler(22*degree,0*degree,0*degree,cs{2},ssO),...
  orientation.goss(cs{2},ssO),...
  orientation.byEuler(70*degree,45*degree,0*degree,cs{2},ssO),...
  orientation.byEuler(59*degree,37*degree,63*degree,cs{2},ssO),...
  };

spread = 10 * degree;

grainsComps = cell(1, 7);
for i=1:length(idealOris)
    ori_name = labels{i};
    ori = idealOris{i};
    grainsComps{i} = findByOrientation(grains('indexed'), ori.symmetrise, spread);
end
% 
% Create cell array in grain2d object
grains.prop.texturecomponent = cell(length(grains),1);

% Assign a texture component to each grain (get rid of overlap!)
cs = grains(2).meanOrientation.CS;
ss = grains(2).meanOrientation.SS;
csal = cs;

defaultOri = orientation.byEuler(1*degree, 0, 0, csal, ssO);
grains.prop.idealOri = repmat(defaultOri, length(grains), 1);

h = waitbar(0, 'Assigning an ideal texture component to each grain');
% Loop over all given ideal texture components

%%% IN CASE I SCREW UP THE CODE
% for i = 1:length(idealOris)
%     waitbar(i/size(idealOris, 2))
% 
%     % Loop over all grains within the current texture component
%     for j = 1:length(grainsComps{i})
%         % Get id of current grain to select the current grain
%         condition = grains.phase==0 & grains.id==grainsComps{i}(j).id;
%         grainJ = grains(condition);
%         id = grainsComps{1,i}(j).id;
% 
%         % Check if grain already have been assigned an ideal
%         % orientation
%         if grainJ.idealOri ~= defaultOri
%             % Get mean orientation of grain and give specimen symmetry
%             grainJOri = orientation(grainJ.meanOrientation,csal,ssO);
%             % Angle between new ideal orientation and its mean
%             % orientation
%             newAngle = angle(idealOris{i},grainJOri);
%             % Angle between current ideal orientation and its mean
%             % orientation
%             oldAngle = angle(grainJ.idealOri,grainJOri);
%             % Give new ideal orientation if new angle is lower than
%             % old angle
%             if newAngle < oldAngle
%                 grains(condition).idealOri = idealOris{i};
%                 %ERLING
%                 grains.prop.texturecomponent{id} = labels{i};
%                 %ERLING
%             end
%         % Give ideal orientation if already not given
%         else
%             grains(condition).idealOri = idealOris{i};
%             %ERLING
%             grains.prop.texturecomponent{id} = labels{i};
%             %ERLING
%         end
%     end
% end
% close(h)
%%% IN CASE I SCREW UP THE CODE

for i = 1:length(idealOris)
    waitbar(i/size(idealOris, 2))

    % Loop over all grains within the current texture component
    for j = 1:length(grainsComps{1,i})
        % Get id of current grain to select the current grain
        condition = grains.phase==0 & grains.id==grainsComps{1,i}(j).id;
        grainJ = grains(condition);
        id = grainsComps{1,i}(j).id;

        % Check if grain already have been assigned an ideal
        % orientation
        if grainJ.idealOri ~= defaultOri
            % Get mean orientation of grain and give specimen symmetry
            grainJOri = orientation(grainJ.meanOrientation,csal,ssO);
            % Angle between new ideal orientation and its mean
            % orientation
            newAngle = angle(idealOris{i},grainJOri);
            % Angle between current ideal orientation and its mean
            % orientation
            oldAngle = angle(grainJ.idealOri,grainJOri);
            % Give new ideal orientation if new angle is lower than
            % old angle
            if newAngle < oldAngle
                grains(condition).idealOri = idealOris{i};
                %ERLING
                grains.prop.texturecomponent{id} = labels{i};
                %ERLING
            end
        % Give ideal orientation if already not given
        else
            grains(condition).idealOri = idealOris{i};
            %ERLING
            grains.prop.texturecomponent{id} = labels{i};
            %ERLING
        end
    end
end
close(h)




% Regroup grains with same ideal orientation
grainsBr = grains(grains.idealOri==idealOris{1});
grainsCu = grains(grains.idealOri==idealOris{2});
grainsCube = grains(grains.idealOri==idealOris{3});
grainscubeND = grains(grains.idealOri==idealOris{4});
grainsGoss = grains(grains.idealOri==idealOris{5});
grainsP = grains(grains.idealOri==idealOris{6});
grainsS = grains(grains.idealOri==idealOris{7});

% Volume fraction from grains' mean orientation
MgBr = 100*sum(grainsBr.area)/sum(grains('al').area);
MgCu = 100*sum(grainsCu.area)/sum(grains('al').area);
MgCube = 100*sum(grainsCube.area)/sum(grains('al').area);
MgcubeND = 100*sum(grainscubeND.area)/sum(grains('al').area);
MgGoss = 100*sum(grainsGoss.area)/sum(grains('al').area);
MgP = 100*sum(grainsP.area)/sum(grains('al').area);
MgS = 100*sum(grainsS.area)/sum(grains('al').area);
%
sorted_grains = {grainsBr,grainsCu,grainsCube,grainscubeND,grainsGoss,grainsP,grainsS};
sorted_grains_strings = {'grainsBr','grainsCu','grainsCube','grainscubeND','grainsGoss','grainsP','grainsS'};
frequencies = [];
for i=1:length(sorted_grains)
    frequencies = [frequencies,length(sorted_grains{i})];
end
volume_labels = {'Br','Cu','Cube','CubeND','Goss','P','S'};
volume_fractions = {label,MgBr,MgCu,MgCube,MgcubeND,MgGoss,MgP,MgS};
vol_exp = [volume_fractions];
nums_exp = [label num2cell(frequencies)];
%% 
writecell(vol_exp,'/Users/erlingaaresvarli/Library/Mobile Documents/com~apple~CloudDocs/Masteroppgave/OneDrive - NTNU/MATLAB/Resultater - MATLAB/Component-fractions/Newcompfractions/Volfrac.xlsx','WriteMode','append');
writecell(nums_exp,'/Users/erlingaaresvarli/Library/Mobile Documents/com~apple~CloudDocs/Masteroppgave/OneDrive - NTNU/MATLAB/Resultater - MATLAB/Component-fractions/Newcompfractions/Frequencies.xlsx','WriteMode','append');
% Saving grains
save(append('/Users/erlingaaresvarli/Library/Mobile Documents/com~apple~CloudDocs/Masteroppgave/OneDrive - NTNU/MATLAB/Resultater - MATLAB/Grain objects with texture components/',label,'-grains'),'grains')

% Plot grains with boundaries

plots = [];
legends = {};
grainsToPlot = {grainsBr, grainsCu, grainsCube, grainscubeND,...
    grainsGoss, grainsP, grainsS, grains('notIndexed')};
grainsToPlotLabels = labels;
grainsToPlotLabels{end + 1} = 'notIndexed';
grainsColours = {'b', 'g', 'r', [1 0.55 0], 'y', 'c', 'm', 'k'};

figure
for i = 1:length(grainsToPlot)
    theseGrains = grainsToPlot{i};
    theseGrainsLabel = grainsToPlotLabels{i};
    theseGrainsColour = grainsColours{i};
    if ~isempty(theseGrains)
        p = plot(theseGrains, theseGrains.area, 'facecolor',...
            theseGrainsColour);
        hold on
        plots = [plots p(1)];
        legends = [legends; {theseGrainsLabel}];
    end
end

gb = grains.boundary('al', 'al');
hold on
plot(gb, 'linecolor', [0 0 0], 'linewidth', 1)
%legends = [legends; {'notIndexed'}];
%legend(plots, legends, 'orientation', 'horizontal');
%legend(grainsToPlotLabels)
legend(plots,legends,'orientation','vertical','location','east')
hold off

%% Grain classification
%print(fullfile('outpath',specimen,filenameom),'-dpdf','-bestfit')
setMTEXpref('figSize','huge')
setMTEXpref('FontSize',30)


%Calculate fraction recrystallized and write to file
%outpath = '/Users/erlingaaresvarli/Library/Mobile Documents/com~apple~CloudDocs/Masteroppgave/MATLAB/FracHab-FracRX/';
resultpath = '/Users/erlingaaresvarli/Library/Mobile Documents/com~apple~CloudDocs/Masteroppgave/OneDrive - NTNU/MATLAB/Resultater - MATLAB/XHab/Back-annealed';

[grains, grainsOR, grainsSub, grainsRex, graindata] = ebsd_fraction_recrystallized_with_og(grains,resultpath,label);
Xrex = sum(grains(grains.RX == 1).area) / sum(grains.area);

write_grains_path = '/Users/erlingaaresvarli/Library/Mobile Documents/com~apple~CloudDocs/Masteroppgave/OneDrive - NTNU/MATLAB/Resultater - MATLAB/Stored grain data/Back-annealed';
ebsd_write_grains_data_modified(grains, fullfile(write_grains_path,append(label,'-grains')),'ebsd',ebsd_clean);
ebsd_write_grains_data_modified(grainsSub, fullfile(write_grains_path,append(label,'-Subgrains')),'ebsd',ebsd_clean);
ebsd_write_grains_data_modified(grainsRex, fullfile(write_grains_path,append(label,'-RXgrains')),'ebsd',ebsd_clean);
ebsd_write_grains_data_modified(grainsOR, fullfile(write_grains_path,append(label,'-ORgrains')),'ebsd',ebsd_clean);
  
% ODF PART
cs = {'notIndexed', crystalSymmetry('m-3m', [4.04 4.04 4.04], 'mineral', 'al')};
ss = specimenSymmetry('1'); % Triclinic
ssO = specimenSymmetry('orthorhombic');

setMTEXpref('FontSize',18)
setMTEXpref('defaultColorMap','white2blackColorMap');
setMTEXpref('markerSize',12)
setMTEXpref('figSize','small')

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
%
% Pole figures
pfAnnotations = @(varargin) text([vector3d.X, vector3d.Y],...
    {'RD', 'TD'}, 'BackgroundColor', 'w', 'tag', 'axesLabels', varargin{:});
setMTEXpref('pfAnnotations', pfAnnotations);
setMTEXpref('defaultColorMap','white2blackColorMap');
%

levelsPDF = [0,1,2,4,6,8];
figure
plotPDF(odf3, h, 'upper', 'projection', 'eangle', 'contourf',levelsPDF)
CLim(gcm,[0 8]);

print(fullfile('/Users/erlingaaresvarli/Library/Mobile Documents/com~apple~CloudDocs/Masteroppgave/OneDrive - NTNU/MATLAB/Resultater - MATLAB/Pole figures/Micro',append(label,'-polefig')),'-dpng');

% Exporting ODF and grain data
odfresultpath = '/Users/erlingaaresvarli/Library/Mobile Documents/com~apple~CloudDocs/Masteroppgave/OneDrive - NTNU/MATLAB/Resultater - MATLAB/ODF data/Micro';
fnameodf = fullfile(odfresultpath, append(label,'-ODFdata.mat'));
save(fnameodf,'odf')

% phi2 sections
setMTEXpref('FontSize',25)
setMTEXpref('figSize','medium')
setMTEXpref('defaultColorMap','white2blackColorMap');
setMTEXpref('markerSize',12)

levelsODF = [0,1,2,4,8,12,15,20,25,30];
%
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

% Calculating area weighed grain sizes
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

Xsub = sum(grains(grains.RX == 0).area) / sum(grains.area);

% Fibres and data export

%Exporting grain data
output_data = {label,average_grain_size_Sub,average_grain_size_RX,average_grain_size_OG,Xrex,Xsub};
writecell(output_data,'/Users/erlingaaresvarli/Library/Mobile Documents/com~apple~CloudDocs/Masteroppgave/OneDrive - NTNU/MATLAB/Resultater - MATLAB/General-results-GOS 2 degrees.xlsx','WriteMode','append');


% Writing fibres to file 
f = fibre(cu, br, cs, ssO);
fibrepath_beta = '/Users/erlingaaresvarli/Library/Mobile Documents/com~apple~CloudDocs/Masteroppgave/OneDrive - NTNU/MATLAB/Resultater - MATLAB/Fibredata/Beta/Micro';

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
fibrepath_cubegoss = '/Users/erlingaaresvarli/Library/Mobile Documents/com~apple~CloudDocs/Masteroppgave/OneDrive - NTNU/MATLAB/Resultater - MATLAB/Fibredata/Cube-Goss/Micro';

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
fibrepath_cubecubeND45 = '/Users/erlingaaresvarli/Library/Mobile Documents/com~apple~CloudDocs/Masteroppgave/OneDrive - NTNU/MATLAB/Resultater - MATLAB/Fibredata/Cube-CubeND45/Micro';

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
