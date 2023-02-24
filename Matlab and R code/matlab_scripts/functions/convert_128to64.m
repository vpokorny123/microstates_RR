function [interpEEG] = convert_128to64(EEG)
addpath('/labs/srslab/static_files/shared_apps/matlab_toolboxes/ssk_eegtoolbox/ICAcleanEEG.v.1.3/')%<<-- change
addpath('/labs/srslab/static_files/shared_apps/matlab_toolboxes/ssk_eegtoolbox/ICAcleanEEG.v.1.3/sphspline0.2/')%<<-- change
addpath('/labs/srslab/static_files/shared_apps/matlab_toolboxes/ssk_eegtoolbox/ICAcleanEEG.v.1.3/montage/')%<<-- change
addpath('/labs/srslab/static_files/shared_apps/matlab_toolboxes/ssk_eegtoolbox/ICAcleanEEG.v.1.3/preprocess/')%<<-- change


BS128 = importdata('/labs/srslab/data_main/Dichotic_NPGGTFPENS/scripts/matfiles/BS128_montage.mat'); %<<-- change
%BS128Ei = BS128.Ei;


BS64 = importdata('/labs/srslab/data_main/Dichotic_NPGGTFPENS/scripts/matfiles/BS64_montage.mat');%<<-- change
BS64Ei = BS64.Ei;
%

nElects = size(EEG.data2interp, 1);
Ei = zeros(nElects, 3);

Ei = BS128.Ei;
data = EEG.data2interp;
%% channel interpolation
stype = 'sp';
lambda = 1e-8;
n = 50;
m = 3;

%SEi - biosemi; VEi - brainvision
% convert biosemi to brainvision
interpEEG = [];
for channel = 1:length(BS64Ei)
    tic
    F = BS64Ei(channel,:); %  Cartesian x,y,z coordinate of the electrodes to interpolate (F)
    [Ginv, g, G] = sserpgfcn(Ei, F, stype, lambda, n, m);
    %     fullElecData = zeros(size(Ei,1), size(data, 2));
    fullElecData = [];
    fullElecData(:,:) = data;
    C = sserpweights(fullElecData, Ginv); % <-------------- Need to use full-electrode data for sserpweights.m [Updated 8/17/2015].
    % C = sserpweights(tecnt_af.data, Ginv);
    interpEEG(channel, :) = sserp(C, g, stype);
    clear Ginv g G C
    disp(['done interpolating channel ', num2str(channel), ' in ' num2str(toc)] )
end

% topoplot_besa_mod4eColor(interpEEG, 'E:\programs\ICAcleanEEG.v1.4beta\montage\BioSemi_128_elecN.elp','maplimits', [min(interpEEG) max(interpEEG)] )