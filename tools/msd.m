% analysis MSD for diffusion validation
clear all
clc
t0 = 0.00226649 %second
l0 = 0.077 %um
% 10dna_21um
table_21um = importdata('output_data_10dna_21um_2.txt');
t_21um = table_21um.data(:,2) * t0;
msd_x_21um=table_21um.data(:,3) * l0^2;
msd_y_21um=table_21um.data(:,4) * l0^2;
msd_xy_21um = msd_x_21um + msd_y_21um;
figure(1)
plot(t_21um,msd_xy_21um,'b')
hold on
%plot([0:0.01:3],1.28*[0:0.01:3],'r')
%fit_21um=polyfit(t_21um(300:end),msd_xy_21um(300:end),1)
f = fittype('4*D*x');
[fit_21um,gof_21um,fitinfo_21um] = fit(t_21um(100:end),msd_xy_21um(100:end),f,'StartPoint',[1]);
plot(t_21um(100:end),fit_21um.D*4*t_21um(100:end),'b')


% 10dna_21um_2
table_21um_2 = importdata('output_21um_10chains_newcode.txt');
t_21um_2 = table_21um_2.data(:,2) * t0;
msd_x_21um_2=table_21um_2.data(:,3) * l0^2;
msd_y_21um_2=table_21um_2.data(:,4) * l0^2;
msd_xy_21um_2 = msd_x_21um_2 + msd_y_21um_2;
plot(t_21um_2,msd_xy_21um_2,'bo')
hold on
plot([0:0.01:3],1.28*[0:0.01:3],'r')


% % 5dna_42um
% 
% table_42um = importdata('output_data_5dna_42um.txt');
% t_42um = table_42um.data(:,2) * t0;
% msd_x_42um=table_42um.data(:,3) * l0^2;
% msd_y_42um=table_42um.data(:,4) * l0^2;
% msd_xy_42um = msd_x_42um + msd_y_42um;
% plot(t_42um,msd_xy_42um,'k')
% hold on
% [fit_42um,gof_42um,fitinfo_42um] = fit(t_42um(100:end),msd_xy_42um(100:end),f,'StartPoint',[1]);
% plot(t_42um(100:end),fit_42um.D*4*t_42um(100:end),'k')
% 
% % % 10dna_42um
% % table_42um_10dna = importdata('output_data_10dna_42um.txt');
% % t_42um_10dna = table_42um_10dna.data(:,2) * t0;
% % msd_x_42um_10dna=table_42um_10dna.data(:,3) * l0^2;
% % msd_y_42um_10dna=table_42um_10dna.data(:,4) * l0^2;
% % msd_xy_42um_10dna = msd_x_42um_10dna + msd_y_42um_10dna;
% % plot(t_42um_10dna,msd_xy_42um_10dna,'k')
% % hold on
% % plot([0:0.01:3],0.60*[0:0.01:3],'r')
% 
% % 5dna_42um_2
% 
% % table_42um_2 = importdata('output_data_5dna_42um_2.txt');
% % t_42um_2 = table_42um_2.data(:,2) * t0;
% % msd_x_42um_2=table_42um_2.data(:,3) * l0^2;
% % msd_y_42um_2=table_42um_2.data(:,4) * l0^2;
% % msd_xy_42um_2 = msd_x_42um_2 + msd_y_42um_2;
% % plot(t_42um_2,msd_xy_42um_2,'k*')
% % %plot(log10(t_42um),log10(msd_xy_42um),'k')
% % hold on
% % plot([0:0.01:3],0.60*[0:0.01:3],'r')
% 
% % 10dna_84um
% table_84um = importdata('output_data_10dna_84um.txt');
% t_84um = table_84um.data(:,2) * t0;
% msd_x_84um=table_84um.data(:,3) * l0^2;
% msd_y_84um=table_84um.data(:,4) * l0^2;
% msd_xy_84um = msd_x_84um + msd_y_84um;
% plot(t_84um,msd_xy_84um,'g')
% hold on
% [fit_84um,gof_84um,fitinfo_84um] = fit(t_84um(100:end),msd_xy_84um(100:end),f,'StartPoint',[1]);
% plot(t_84um(100:end),fit_84um.D*4*t_84um(100:end),'g')
% hold off
% % 
% % % 
% % % % 10dna_84um_2
% % % table_84um_2 = importdata('output_data_10dna_84um_2.txt');
% % % t_84um_2 = table_84um_2.data(:,2) * t0;
% % % msd_x_84um_2=table_84um_2.data(:,3) * l0^2;
% % % msd_y_84um_2=table_84um_2.data(:,4) * l0^2;
% % % msd_xy_84um_2 = msd_x_84um_2 + msd_y_84um_2;
% % % %plot(t_84um_2,msd_xy_84um_2,'g*')
% % % hold on
% % % %plot([0:0.01:3],0.60*[0:0.01:3],'r')
% % 
% % 
% % 
% % 
% figure(2)
% Rg=[0.594,0.85,1.43];
% H=2;
% Rg_H=Rg/H;
% D=[fit_21um.D,fit_42um.D,fit_84um.D]
% D2=[1.28,0.72,0.3]/4
% Dbulk=[0.53,0.37,0.22];
% D_Dbulk=D./Dbulk
% plot((Rg_H),(D_Dbulk),'b*')
% hold on
% fit_scale_type=fittype('a*x')
% [fit_scale,gof_scale,fitinfo_scale] = fit(Rg_H',D_Dbulk',fit_scale_type,'StartPoint',[1]);
% plot((Rg_H(1):0.01:Rg_H(3)),(0.3*(Rg_H(1):0.01:Rg_H(3)).^(-0.7)),'r')
% hold off
