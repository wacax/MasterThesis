function str  = EEGLoad(Participant, VP)

global BBCI_PRINTER
BBCI_PRINTER = 1;

figdir = 'd:\Wacax\TU Berlin\figures\'
% 
% amplitudesCz = zeros(1, 3);
% amplitudesTP7 = zeros(1, 3);
% amplitudesTP8= zeros(1, 3);
% latencyCz = zeros(1, 3);
% latencyTP7 = zeros(1, 3);
% latencyTP8 = zeros(1, 3);
% 
% numbers = {'' '02' '03'}

erp_av = {};
erp_r_av = {};

for i = 1:length(Participant),
%scaled by 100
file = sprintf('%s%s%s%s', Participant{i}, '\navigation_train_std_', VP{i},'*');
hdr= eegfile_readBVheader(file);
Wps= [40 49]/hdr.fs*2;
[n, Ws]= cheb2ord(Wps(1), Wps(2), 3, 50);
[filt.b, filt.a]= cheby2(n, 50, Ws);
[cnt, mrk] = eegfile_loadBV(file, 'fs', 100, 'filt', filt, 'filtType', 2);

% Define some settings
disp_ival= [-200 1000];
ref_ival= [-200 0];
crit_maxmin= 70;
crit_ival= [100 800];
crit_clab= {'F9,z,10','AF3,4'};
clab= {'Cz','PO7'};
colOrder= [1 0 1; 0.4 0.4 0.4];

% high-pass filtering to reduce drifts
b= procutil_firlsFilter(0.5, cnt.fs);
cnt= proc_filtfilt(cnt, b);

%Define Important Markers
mrkDef = {[11:19],[1:9];'Target','Non-target'}%Which ones are the targets and which ones are not

%Processed Marker Structure
mrk = mrk_defineClasses(mrk, mrkDef)

% Artifact rejection based on variance criterion
mrk= reject_varEventsAndChannels(cnt, mrk, disp_ival, 'verbose', 1);

% Here, the first class corresponds to the presentation of target stimuli. You get in indices of target markers by
it= find(mrk.y(1,:));
it(1:10)
% The following displays the scalp distribution at the time point of the first target presentation:
%scalpPlot(mnt, cnt.x(mrk.pos(it(1)),:)) % This function is deprecated
%apparently

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Naw they be workin'

% segmentation of continuous data in 'epochs' based on markers
epo= cntToEpo(cnt, mrk, disp_ival); %ask Martijn which one is better
epo= proc_rejectArtifactsMaxMin(epo, crit_maxmin, 'clab',crit_clab,...
'ival',crit_ival, 'verbose',1); %this one is still working

%Most Important Numbers
% Baseline subtraction, and calculation of a measure of discriminability
epo= proc_baseline(epo, ref_ival); %still working
% epo_r= proc_r_square_signed(epo); %Still!!!
epo_r = proc_rocAreaValues(epo);

%Creating the montage based on epochs' channels
mnt= getElectrodePositions(epo.clab);
mnt = mnt_adaptMontage(mnt, epo.clab);
%TO DO inlude the montage preprocessing

% Select some discriminative intervals, with constraints to find N2, P2, P3 like components.
fig_set(1);
constraint= ...
       {{-1, [100 300], {'I#','O#','PO7,8','P9,10'}, [50 300]}, ...
       {1, [200 350], {'P3-4','CP3-4','C3-4'}, [200 400]}, ...
       {1, [400 500], {'P3-4','CP3-4','C3-4'}, [350 600]}};
[ival_scalps, nfo]= ...
     select_time_intervals(epo_r, 'visualize', 1, 'visu_scalps', 1, ...
                           'title', untex(file), ...
                           'clab',{'not','E*'}, ...
                           'constraint', constraint); %it works fuck yea
                       
close 
                       
%printFigure('r_matrix', [18 13]); %This command from the wiki doesn't work
ival_scalps= visutil_correctIvalsForDisplay(ival_scalps, 'fs',epo.fs);
printFigure([figdir 'ivals' VP{i}], [19 12], 'format', 'svg');

fig_set(3)
H= grid_plot(epo, mnt, defopt_erps, 'colorOrder',colOrder);
%grid_addBars(epo_r, 'h_scale',H.scale);
%IMPORTANT adding bars is not possible due to the lack of field scale in
%structure H
printFigure([figdir 'erp' VP{i}], [19 12]);

close 

fig_set(2);
H= scalpEvolutionPlusChannelPlusRsquared(epo,epo_r, mnt, clab, ival_scalps, defopt_scalp_erp2, ...
                             'colorOrder',colOrder);
grid_addBars(epo_r);
printFigure([figdir 'erp_topo' VP{i}], [20  4+5*size(epo.y,1)]);

close 
% fig_set(4, 'shrink',[1 2/3]);
% scalpEvolutionPlusChannel(epo_r, mnt, clab, ival_scalps, defopt_scalp_r2);
% printFigure(['erp_topo_r'], [20 9]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%OR JUST SIMPLY USE JOHANNES' PLOTTING FUNCTION, it does all heavy lifting
%[arf1 arf2 arf3] = stdERPplots(epo, epo_r)
%input all the extra arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5

% collect averages for calculating grand averages later
%Not ready yet
erp_av = sprintf('%s%i', 'erp_av', i);
erp_r_av = sprintf('%s%i', 'erp_r_av', i);

str.(erp_av) = proc_average(epo);
str.(erp_r_av) = epo_r;
%   
% Extract and Store ERP components
% derp = erp_components(epo, 'p', [200 300], 'Cz');
% amplitudesCz(i) = derp.amplitude;
% latencyCz(i) = derp.latency;
% derp2 = erp_components(epo, 'p', [200 300], 'TP7');
% amplitudesTP7(i) = derp2.amplitude;
% latencyTP7(i) = derp2.latency;
% derp3 = erp_components(epo, 'p', [200 300], 'TP8');
% amplitudesTP8(i) = derp3.amplitude;
% latencyTP8(i) = derp3.latency;

end
% 
% str.dataav = erp_av
% str.datarav = erp_r_av
% 
% erp_av = proc_grandAverage(erp_av);
% erp_r_av = proc_grandAverage(erp_r_av);
% 
% str.grandav = erp_av
% str.grandrav = erp_r_av
% 
% 
% stdERPplots(erp_av, erp_r_av)
% 
% amplitudesCz = mean(amplitudesCz)
% latencyCz = mean(latencyCz)
% amplitudesTP7 = mean(amplitudesTP7)
% latencyTP7 = mean(latencyTP7)
% amplitudesTP8 = mean(amplitudesTP8)
% latencyTP8 = mean(latencyTP8)

end