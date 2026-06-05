%% This code infers temperatures from seismic velocities

clear all;
close all;

%%---------------
%% USER INPUTS needed: composition (comp) where:
% compositions are Ha- harzburgite, Py-pyrolite, Ol-olivine, Ba-basalt
comp = 'Ha';


%%---------------
%% DIRECTORIES needed where the HeFESTo data are store, and where to save the extracted quantities from HeFESTo.
%% THESE have to be changed.

directory = 'D:\HPDesktop\PhD\Research\ForPublications\Paper1\phases\';
directory1 = 'D:\HPDesktop\PhD\Research\ForPublications\Paper1\lookup2011\Interp_1K_interval\';

%%for saving the extracted quantities
directory0 = 'D:\HPDesktop\PhD\Research\ForPublications\Paper3\Results\';

%%-----------------------
%% FILENAMES of the HeFESTo generated data that contains calculated velocities at certain P and T for specific composition.
%% The user can use the following data (.56) instead of the interpolated lookup data (with a filename
%%InterpolatedVelocities_allP&T)

fname = [{'olivine.56'},{'harzburgite.56'},{'pyrolite.56'}, {'basalt.56'}];


if comp == 'Ol'
    fname = {fname{1}};
elseif comp == 'Ha'
    fname = {fname{2}};
elseif comp == 'Py'
    fname = {fname{3}};
elseif comp == 'Ba'
    fname = {fname{4}};
end



%% --------------------------------------------------------------------------------
%% Here the USER needs to load the data VELOCITY PROFILE with depth and velocity values.
%% Note that this has to be changed but preferrably maintain the naming, which will be used in extracting quatities from HeFESTo %% data.
%% This part is to extract from the actual data the velocity ranges (low and high velocities) and their corresponding depths

%%% FOR cratonic domains REGIONAL seismic data (Norway-Sweden-Finland)
%source = 'Maupin2022_Bruneton2004_ani';
%dir = 'D:\HPDesktop\PhD\Research\ForPublications\Paper3\InputData\';
%Name = strcat('Vsv_regionalData_',source,'.mat');
%data = load(strcat(dir,Name));

% region = 'NS'; %Norway-Sweden = NS, Finland = F, cratonic domains = CD

%% for average and velocity range
%%if region == 'CD'
%    %for OVERALL CRATONIC DOMAINS
%    VsR_T = data.Fennoscandia(:,2);
%    VsR_Fi = data.Fennoscandia(:,5);
%    Vs_lo = 0.75.*VsR_T + 0.25.*VsR_Fi;
%    Vs_up = Vs_lo;
%    z_data = data.Fennoscandia(:,1);

%    elseif region == 'NS'
%    %for Sweden-Norway (Maupin et al. 2022)
%    Vs_lo = data.Vsv(:,7);%data.Vsv(:,5); column 5 - Vsv; col 6 - Vsh; col 7 - Vs(voigt)
%    Vs_up = data.Vsv(:,7);
%    z_data = data.Vsv(:,1);

%    elseif region == 'F'
%    %for Finland (Bruneton et al. 2004)
%    Vs_lo = data.Fennoscandia(:,5);
%    Vs_up = data.Fennoscandia(:,5);
%    z_data = data.Fennoscandia(:,1);
%end

%Vs_ave = 0.5.*(Vs_up+Vs_lo);
%Vs_err = 0.5.*(Vs_up-Vs_lo);




%% ----------------------------
%% EXTRACTING PROPERTIES for specific pressure or depth from interpolated temperatures data
% The loaded data, vec, can be changed to uninterpolated data like harzburgite.56 (see how it looks like)

vec = load(strcat(directory1,'InterpolatedVelocities_allP&T_',string(fname),'_Qcor.mat'));

%% reference T, z, P
T_ = vec.T_int';
z_ = vec.z;
P_ = vec.P;

if comp == 'Ha' | comp == 'Py'
    vec_p = load(strcat(directory,'InterpolatedPhases_allP&T_',string(fname),'.mat'));
else
    vec_p = [];
end


%% A. SEISMICALLY INFERRED TEMPERATURES
%% NOTE: The z_data is the depth values from the seimic velocity profile loaded by the user. So check the path

%% Getting the nearest z value and its index from the HeFESTo to represent the user's data
z_rep = repmat(z_, [1 length(z_data)]);
[minVal, closestIndex] = min(abs(z_rep-z_data'));
closestVal_z = z_(closestIndex);
diff = z_data - closestVal_z;
Indx_fi = find(abs(diff)<=1.5);
closestVal_z_fi = closestVal_z(Indx_fi);
closestIndex_fi = closestIndex(Indx_fi);

%% EXTRACTING average velocity and velocity ranges

for k=1:length(closestIndex_fi)

    Ind = closestIndex_fi(k);
    P(k,1) = P_(Ind);
    
    [diff_av(k), Ind_Vs_ave(k)] = min(abs(round((vec.VsQ_int(Ind,:)-Vs_ave(Indx_fi(k))),4)));    
    [diff_lo(k), Ind_Vs_lo(k)] = min(abs(round((vec.VsQ_int(Ind,:)-Vs_lo(Indx_fi(k))),4)));       
    [diff_up(k), Ind_Vs_up(k)] = min(abs(round((vec.VsQ_int(Ind,:)-Vs_up(Indx_fi(k))),4)));   
    
    T1(k,1) = T_(Ind_Vs_up(k));
    T2(k,1) = T_(Ind_Vs_lo(k));
    Tav(k,1) = T_(Ind_Vs_ave(k));
    
%%----Vs calculated from HeFESTo that is closest to the inputted data---

    Vs_ave_HeF(k,1) = vec.VsQ_int(Ind,Ind_Vs_ave(k)); 
    Vs_up_HeF(k,1) = vec.VsQ_int(Ind,Ind_Vs_up(k));
    Vs_lo_HeF(k,1) = vec.VsQ_int(Ind,Ind_Vs_lo(k));
end




%% B. EXTRACTING PHASE PROPORTIONS

if comp == 'Ha' | comp == 'Py'

    %% Method 1: averaging the phase proportions to represent for average temperatures

    for k=1:length(closestIndex_fi)

        Ind = closestIndex_fi(k);
        P(k,1) = P_(Ind);

        plg_T1(k) = vec_p.plg_in(Ind,Ind_Vs_up(k));
        plg_T2(k) = vec_p.plg_in(Ind,Ind_Vs_lo(k));
        plg_ave(k) = 0.5*(plg_T1(k) + plg_T2(k));

        sp_T1(k) = vec_p.sp_in(Ind,Ind_Vs_up(k));
        sp_T2(k) = vec_p.sp_in(Ind,Ind_Vs_lo(k));
        sp_ave(k) = 0.5*(sp_T1(k) + sp_T2(k));

        opx_T1(k) = vec_p.opx_in(Ind,Ind_Vs_up(k));
        opx_T2(k) = vec_p.opx_in(Ind,Ind_Vs_lo(k));
        opx_ave(k) = 0.5*(opx_T1(k) + opx_T2(k));

        c2c_T1(k) = vec_p.c2c_in(Ind,Ind_Vs_up(k));
        c2c_T2(k) = vec_p.c2c_in(Ind,Ind_Vs_lo(k));
        c2c_ave(k) = 0.5*(c2c_T1(k) + c2c_T2(k));

        cpx_T1(k) = vec_p.cpx_in(Ind,Ind_Vs_up(k));
        cpx_T2(k) = vec_p.cpx_in(Ind,Ind_Vs_lo(k));
        cpx_ave(k) = 0.5*(cpx_T1(k) + cpx_T2(k));

        gt_T1(k) = vec_p.gt_in(Ind,Ind_Vs_up(k));
        gt_T2(k) = vec_p.gt_in(Ind,Ind_Vs_lo(k));
        gt_ave(k) = 0.5*(gt_T1(k) + gt_T2(k));

        cpv_T1(k) = vec_p.cpv_in(Ind,Ind_Vs_up(k));
        cpv_T2(k) = vec_p.cpv_in(Ind,Ind_Vs_lo(k));
        cpv_ave(k) = 0.5*(cpv_T1(k) + cpv_T2(k));

        ol_T1(k) = vec_p.ol_in(Ind,Ind_Vs_up(k));
        ol_T2(k) = vec_p.ol_in(Ind,Ind_Vs_lo(k));
        ol_ave(k) = 0.5*(ol_T1(k) + ol_T2(k));

        wa_T1(k) = vec_p.wa_in(Ind,Ind_Vs_up(k));
        wa_T2(k) = vec_p.wa_in(Ind,Ind_Vs_lo(k));
        wa_ave(k) = 0.5*(wa_T1(k) + wa_T2(k));

        ri_T1(k) = vec_p.ri_in(Ind,Ind_Vs_up(k));
        ri_T2(k) = vec_p.ri_in(Ind,Ind_Vs_lo(k));
        ri_ave(k) = 0.5*(ri_T1(k) + ri_T2(k));

        qtz_T1(k) = vec_p.qtz_in(Ind,Ind_Vs_up(k));
        qtz_T2(k) = vec_p.qtz_in(Ind,Ind_Vs_lo(k));
        qtz_ave(k) = 0.5*(qtz_T1(k) + qtz_T2(k));

        ky_T1(k) = vec_p.ky_in(Ind,Ind_Vs_up(k));
        ky_T2(k) = vec_p.ky_in(Ind,Ind_Vs_lo(k));
        ky_ave(k) = 0.5*(ky_T1(k) + ky_T2(k));
    end

    phasePercent1 = [ol_ave',opx_ave',cpx_ave',gt_ave',plg_ave',sp_ave',c2c_ave',cpv_ave',qtz_ave',ky_ave',wa_ave',ri_ave'];
    phasePercent_T1 = [ol_T1',opx_T1',cpx_T1',gt_T1',plg_T1',sp_T1',c2c_T1',cpv_T1',qtz_T1',ky_T1',wa_T1',ri_T1'];
    phasePercent_T2 = [ol_T2',opx_T2',cpx_T2',gt_T2',plg_T2',sp_T2',c2c_T2',cpv_T2',qtz_T2',ky_T2',wa_T2',ri_T2'];


    %% Method 2: using the average T

    T_rep = repmat(T_, [1 length(T1)]); %[1 length(z_data)]
    [minValT, closestIndexT] = min(abs(T_rep-Tav'));  % from T_ave
    closestVal_T = T_(closestIndexT);

    for k=1:length(closestIndexT)

        Ind = closestIndex(k);

        plg_ave2(k) = vec_p.plg_in(Ind,closestIndexT(k));
        sp_ave2(k) = vec_p.sp_in(Ind,closestIndexT(k));
        opx_ave2(k) = vec_p.opx_in(Ind,closestIndexT(k));
        c2c_ave2(k) = vec_p.c2c_in(Ind,closestIndexT(k));
        cpx_ave2(k) = vec_p.cpx_in(Ind,closestIndexT(k));
        gt_ave2(k) = vec_p.gt_in(Ind,closestIndexT(k));
        cpv_ave2(k) = vec_p.cpv_in(Ind,closestIndexT(k));
        ol_ave2(k) = vec_p.ol_in(Ind,closestIndexT(k));
        wa_ave2(k) = vec_p.wa_in(Ind,closestIndexT(k));
        ri_ave2(k) = vec_p.ri_in(Ind,closestIndexT(k));
        qtz_ave2(k) = vec_p.qtz_in(Ind,closestIndexT(k));
        ky_ave2(k) = vec_p.ky_in(Ind,closestIndexT(k));
    end

    %% better to use METHOD 2

    phasePercent2 = [ol_ave2',opx_ave2',cpx_ave2',gt_ave2',plg_ave2',sp_ave2',c2c_ave2',cpv_ave2',qtz_ave2',ky_ave2',wa_ave2',ri_ave2'];

else
    phasePercent1 = [ones(length(closestIndex_fi),1),zeros(length(closestIndex_fi),11)];
    phasePercent_T1 = [ones(length(closestIndex_fi),1),zeros(length(closestIndex_fi),11)];
    phasePercent_T2 = [ones(length(closestIndex_fi),1),zeros(length(closestIndex_fi),11)];
    phasePercent2 = [ones(length(closestIndex_fi),1),zeros(length(closestIndex_fi),11)];
end

%% SAVING THE INFERRED QUANTITIES
newMatRange = [z_data(Indx_fi,1), closestVal_z_fi, P, Vs_lo(Indx_fi,1),  Vs_up(Indx_fi,1), T1, T2];
newMatAve = [z_data(Indx_fi,1), P, Vs_ave(Indx_fi,1), Tav];
newMatHeF = [z_data(Indx_fi,1), P, Vs_ave_HeF, Vs_lo_HeF, Vs_up_HeF];

save(strcat(directory0,'SeismicInferences_',comp,'_Qcor_','source,'.mat'),'newMatRange','newMatAve','newMatHeF','phasePercent1','phasePercent2', 'phasePercent_T1', 'phasePercent_T2');




















%%---BELOW IS SAMPLE EXTRACTION using Schaeffer and Lebedev 2013 data----
% %% C. Seismically inferred temperatures for every lat-lon
% 
% %% for Schaeffer & Lebedev 2013 data
% z_rep = repmat(z_, [1 length(z_data)]);
% [minVal, closestIndex] = min(abs(z_rep-z_data'));
% closestVal_z = z_(closestIndex);
% 
% % for velocity at certain lat and lon
% [row col] = size(Vs_latlon);
% for kk=1:col
%     Vs = Vs_latlon(Indx_fi,kk);
%     
% for k=1:length(closestIndex_fi)
%     Ind = closestIndex_fi(k);
% %     P(k,1) = P_(Ind);
%     if kk==1
%     P(k,1) = P_(Ind);
%     end
%     
%     [diff(k,kk), Ind_Vs(k,kk)] = min(abs(round((vec.VsQ_int(Ind,:)-Vs(k)),4)));    
%   
% %     Ind_Vs_ = find(round((vec.VsQ_int(Ind,:)-Vs(k)),4)>=0.0009);%find(abs(vec.Mat_int(:,2)-Vs_lo)<=0.009);
% %     Ind_Vs(k,kk) = Ind_Vs_(end);
%     
%     T(k,kk) = T_(Ind_Vs(k,kk));
%     
%      %Vs calculated from HeFESTo that is closest to the inputted data
%     Vs_HeF(k,kk) = vec.VsQ_int(Ind,Ind_Vs(k,kk)); 
% 
% end
% end
% % 
% % Ind32 = find(z_data==32);
% % Vs_data = [z_data(1:Ind32,1), P, Vs_latlon(1:Ind32,:)];
% % T_inferred = [z_data(1:Ind32,1), P, T];
% % Vs_HeF = [z_data(1:Ind32,1), P, Vs_HeF(1:Ind32,:)];
% Vs_data = [z_data(Indx_fi,1), P, Vs_latlon(Indx_fi,:)];
% T_inferred = [z_data(Indx_fi,1), P, T];
% Vs_HeF = [z_data(Indx_fi,1), P, Vs_HeF(:,:)];
% save(strcat(directory0,'SeismicInferences_LatLon_',geotherm,comp,'_Qcor_',source,'.mat'),'Vs_data','T_inferred','Vs_HeF','latlon');

