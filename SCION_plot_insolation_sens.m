%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% SCION - Spatial Continuous Integration %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% Earth Evolution Model %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Coded by BJW Mills and Yinggang Zhang
%%%% b.mills@leeds.ac.uk  and ygzhang@nigpas.ac.cn
%%%%
%%%% plot sensitivity analysis

%%%%%% define colours
c_mean = [255 132 34]./255 ;
c_std = [255 225 192]./255 ;
c_range = [255 225 192]./255 ;

%%%% Proxy color chart
pc1 = [65 195 199]./255 ;
pc2 = [73 167 187]./255 ;
pc3 = [82 144 170]./255 ;
pc4 = [88 119 149]./255 ;
pc5 = [89 96 125]./255 ;
pc6 = [82 56 100]./255 ;

%%%% output to screen
fprintf('running sens plotting script... \t')
tic

%%%% load geochem data
load('data/geochem_data_2020.mat')

%%%%%%% make figure˙
figure('Color',[0.80 0.80 0.70])

%%%% load geochem data
c = xlsread('cmodel.xlsx','','','basic') ;
c_time = c(:,1)*-1;
c_data = c(:,2);
c_untrend_time = c(:,3);
c_untrend_data = c(:,4);

s = xlsread('smodel.xlsx','','','basic') ;
s_time = s(:,1)*-1;
s_data = s(:,2);
s_untrend_time = s(:,3);
s_untrend_data = s(:,4);

%%%% make column vector
sens.time_myr = sens.time_myr(:,1) ;

plotrange = [-526 -512];


subplot(5,2,1)
hold on
box on
xlabel('Time (Ma)')
ylabel('Relative Size')
%%%% plot this model
plot((sens.time_myr),mean(sens.DEGASS,2),'linewidth',1,'color',c_mean)
plot((sens.time_myr),max(sens.DEGASS,[],2),'linewidth',0.1,'color',c_range)
plot((sens.time_myr),min(sens.DEGASS,[],2),'linewidth',0.1,'color',c_range)

%%%% plot this model
plot((sens.time_myr),mean(sens.BAS_AREA,2),'linewidth',1,'color',c_std)
plot((sens.time_myr),max(sens.BAS_AREA,[],2),'linewidth',0.1,'color',c_range)
plot((sens.time_myr),min(sens.BAS_AREA,[],2),'linewidth',0.1,'color',c_range)

%%%% plot this model
plot((sens.time_myr),mean(sens.GRAN_AREA,2),'linewidth',1,'color',c_std)
plot((sens.time_myr),max(sens.GRAN_AREA,[],2),'linewidth',0.1,'color',c_range)
plot((sens.time_myr),min(sens.GRAN_AREA,[],2),'linewidth',0.1,'color',c_range)

xlim(plotrange)

load('data/Scotese_GAT_2021.mat')
%%%% Temperature
subplot(5,2,2)
hold on
box on
xlabel('Time (Ma)')
ylabel('Temperature')
plot(Scotese_2021_age , Scotese_2021_GAT, 'linewidth',1,'color',pc2)
%%%% plot this model
plot((sens.time_myr),mean(sens.T_gast,2),'linewidth',0.5,'color',c_mean)
plot((sens.time_myr),max(sens.T_gast,[],2),'linewidth',1,'color',c_range)
plot((sens.time_myr),min(sens.T_gast,[],2),'linewidth',1,'color',c_range)


xlim(plotrange)
ylim([0 28])

% phos weathering
subplot(5,2,3)
hold on
box on
xlabel('Time (Ma)')
ylabel('Phosphorus weathering ')
%%%% plot this model
plot((sens.time_myr),mean(sens.phosw,2),'linewidth',0.5,'color',c_mean)
plot((sens.time_myr),max(sens.phosw,[],2),'linewidth',1,'color',c_range)
plot((sens.time_myr),min(sens.phosw,[],2),'linewidth',1,'color',c_range)
xlim(plotrange)

Udata = xlsread('Umodel.xlsx','','','basic') ;
Udata_time = Udata (:,1)*-1;
Udata_data = Udata(:,2);
% anoxia extent in shelf and open ocean
subplot(5,2,4)
hold on
box on
xlabel('Time (Ma)')
ylabel('Relative size')
% %%%% plot this model

yyaxis left
%%%% plot this model
plot((sens.time_myr),mean(sens.ANOX_s,2),'-','linewidth',1,'color',c_mean)
plot((sens.time_myr),max(sens.ANOX_s,[],2),'-','linewidth',0.1,'color',c_range)
plot((sens.time_myr),min(sens.ANOX_s,[],2),'-','linewidth',0.1,'color',c_range)

plot((sens.time_myr),mean(sens.ANOX_i,2),'-','linewidth',1,'color','magenta')
plot((sens.time_myr),max(sens.ANOX_i,[],2),'-','linewidth',0.1,'color','magenta')
plot((sens.time_myr),min(sens.ANOX_i,[],2),'-', 'linewidth',0.1,'color','magenta')


yyaxis right
plot(Udata_time , Udata_data, 'linewidth',1,'color',pc2)
ylim([-1.8 0.1])
xlim(plotrange)

% Marine pyrite S burial in shelf and open ocean
subplot(5,2,5)
hold on
box on
xlabel('Time (Ma)')
ylabel('Marine Organic Carbon Burial (mol/yr)')
%%%% plot this model
plot((sens.time_myr),mean(sens.mocb_i,2),'linewidth',0.5,'color','magenta')
plot((sens.time_myr),max(sens.mocb_i,[],2),'linewidth',1,'color','magenta')
plot((sens.time_myr),min(sens.mocb_i,[],2),'linewidth',1,'color','magenta')

%%%% plot this model
plot((sens.time_myr),mean(sens.mocb_s,2),'linewidth',0.5,'color',c_mean)
plot((sens.time_myr),max(sens.mocb_s,[],2),'linewidth',1,'color',c_range)
plot((sens.time_myr),min(sens.mocb_s,[],2),'linewidth',1,'color',c_range)
xlim(plotrange)

% Marine pyrite S burial in shelf and open ocean
subplot(5,2,6)
hold on
box on
xlabel('Time (Ma)')
ylabel('Marine Pyrite S Burial (mol/yr)')
%%%% plot this model
plot((sens.time_myr),mean(sens.mpsb_i,2),'linewidth',1,'color','magenta')
plot((sens.time_myr),max(sens.mpsb_i,[],2),'linewidth',0.1,'color','magenta')
plot((sens.time_myr),min(sens.mpsb_i,[],2),'linewidth',0.1,'color','magenta')

%%%% plot this model
plot((sens.time_myr),mean(sens.mpsb_s,2),'linewidth',0.5,'color',c_mean)
plot((sens.time_myr),max(sens.mpsb_s,[],2),'linewidth',1,'color',c_range)
plot((sens.time_myr),min(sens.mpsb_s,[],2),'linewidth',1,'color',c_range)
xlim(plotrange)

% Shelf P reservoir size
subplot(5,2,7)
hold on
box on
xlabel('Time (Ma)')
ylabel('Shelf P reservoir size')
%%%% plot this model
plot((sens.time_myr),mean(sens.P_s,2),'linewidth',1,'color',c_mean)
plot((sens.time_myr),max(sens.P_s,[],2),'linewidth',0.1,'color',c_range)
plot((sens.time_myr),min(sens.P_s,[],2),'linewidth',0.1,'color',c_range)
% 
plot((sens.time_myr),mean(sens.P_i,2),'linewidth',1,'color','magenta')
plot((sens.time_myr),max(sens.P_i,[],2),'linewidth',0.1,'color','magenta')
plot((sens.time_myr),min(sens.P_i,[],2),'linewidth',0.1,'color','magenta')
xlim(plotrange)

%%%% d13C record
subplot(5,2,8)
hold on
box on
xlabel('Time (Ma)')
ylabel('\delta^{13}C_{carb}')
%%%% plot data comparison
plot(c_time,c_data,'.','color',pc2)
plot(c_untrend_time,c_untrend_data,'linewidth',0.5,'color',pc2)

%%%% plot this model
plot((sens.time_myr),mean(sens.delta_mccb,2),'linewidth',1,'color',c_mean)
plot((sens.time_myr),max(sens.delta_mccb,[],2),'linewidth',0.1,'color',c_range)
plot((sens.time_myr),min(sens.delta_mccb,[],2),'linewidth',0.1,'color',c_range)
xlim(plotrange)

%%%% d34S record
subplot(5,2,9)
hold on
box on
xlabel('Time (Ma)')
ylabel('\delta^{34}S_{sw}')
%%%% plot data comparison
% plot(s_time,s_data,'.','color',pc2)
plot(s_untrend_time,s_untrend_data,'linewidth',1,'color',pc2)
%%%% plot this model
plot((sens.time_myr),mean(sens.d34s_S,2),'linewidth',1,'color',c_mean)
plot((sens.time_myr),max(sens.d34s_S,[],2),'linewidth',0.1,'color',c_range)
plot((sens.time_myr),min(sens.d34s_S,[],2),'linewidth',0.1,'color',c_range)
xlim(plotrange)

%%%% O2 (%) 
subplot(5,2,10)
hold on
box on
xlabel('Time (Ma)')
ylabel('Atmospheric O_{2} (%)')
%%%% plot this model
plot((sens.time_myr),mean(sens.mrO2.*100,2),'linewidth',1,'color',c_mean)
plot((sens.time_myr),max(sens.mrO2.*100,[],2),'linewidth',0.1,'color',c_range)
plot((sens.time_myr),min(sens.mrO2.*100,[],2),'linewidth',0.1,'color',c_range)

%%%% plot data comparison
for u = 1:2:length(O2_x) - 1
   plot( [O2_x(u) O2_x(u)] , [O2_y(u) O2_y(u+1)] , 'color' , pc2  ) ;     
end

xlim(plotrange)



figure
%%%% CO2ppm
subplot(5,2,1)
set(gca, 'YScale', 'log')
hold on
box on
ylim([100 100000])
xlabel('Time (Ma)')
ylabel('Atmospheric CO_{2} (ppm)')
%%%% plot data comparison
%%%% paleosol
% errorbar(paleosol_age,paleosol_co2,paleosol_low,paleosol_high,'color',[0.4 0.7 0.7],'linestyle','none')
plot(paleosol_age, paleosol_co2,'.','markerfacecolor',pc1,'markeredgecolor',pc1)
%%%% alkenone
% errorbar(alkenone_age,alkenone_co2,alkenone_low,alkenone_high,'color',[0.4 0.7 0.4],'linestyle','none')
plot(alkenone_age, alkenone_co2,'.','markerfacecolor',pc2,'markeredgecolor',pc2)
%%%% boron
% errorbar(boron_age,boron_co2,boron_low,boron_high,'color',[0.7 0.4 0.4],'linestyle','none')
plot(boron_age, boron_co2,'.','markerfacecolor',pc3,'markeredgecolor',pc3)
%%%% stomata
% errorbar(stomata_age,stomata_co2,stomata_low,stomata_high,'color',[0.7 0.7 0.4],'linestyle','none')
plot(stomata_age, stomata_co2,'.','markerfacecolor',pc4,'markeredgecolor',pc4)
%%%% liverwort
% errorbar(liverwort_age,liverwort_co2,liverwort_low,liverwort_high,'color',[0.7 0.7 0.4],'linestyle','none')
plot(liverwort_age, liverwort_co2,'.','markerfacecolor',pc5,'markeredgecolor',pc5)
%%%% phytane
% errorbar(phytane_age,phytane_co2,phytane_low,phytane_high,'color',[0.7 0.7 0.4],'linestyle','none')
plot(phytane_age, phytane_co2,'.','markerfacecolor',pc6,'markeredgecolor',pc6)
%%%% plot this model
plot((sens.time_myr),mean(sens.CO2ppm,2),'linewidth',1,'color',c_mean)
plot((sens.time_myr),max(sens.CO2ppm,[],2),'linewidth',0.1,'color',c_range)
plot((sens.time_myr),min(sens.CO2ppm,[],2),'linewidth',0.1,'color',c_range)
xlim(plotrange)
% 
% %%%% TEMP
subplot(5,2,2)
hold on
box on
xlabel('Time (Ma)')
ylabel('SO42-')
%%%% plot this model
plot((sens.time_myr),mean(sens.SmM,2),'linewidth',1,'color',c_mean)
plot((sens.time_myr),max(sens.SmM,[],2),'linewidth',0.1,'color',c_range)
plot((sens.time_myr),min(sens.SmM,[],2),'linewidth',0.1,'color',c_range)

%%%% plot algeo data window comparison
plot(sconc_max_x,sconc_max_y,'color',pc1)
plot(sconc_min_x,sconc_min_y,'color',pc1)
plot(sconc_mid_x,sconc_mid_y,'color',pc2)
%%%% plot fluid inclusion data comparison
for u = 1:2:length(SO4_x-1)
   plot( [SO4_x(u) SO4_x(u)] , [SO4_y(u) SO4_y(u+1)], 'color' , pc3 ) ;     
end

xlim(plotrange)


%%%%% plotting script finished
fprintf('Done: ')
endtime = toc ;
fprintf('time (s): %d \n', endtime )



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%   Cleanup workspace   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% clear stepnumber
% clear u
% clear numfields
% clear trecords
% clear finalrecord
% clear field_names
% clear n
% clear veclength
% clear xvec
% clear yvec
% clear endtime


%%%%% plotting script finished
fprintf('Done: ')
endtime = toc ;
fprintf('time (s): %d \n', endtime )
