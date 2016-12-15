% function plot_Regimes
%
% LJ YIEW
% Created on  Feb 2016
% Last edited Dec 2016
%
% Plots rafting, nonrafting, no collision regimes as a function of
% wavelength and wave amplitude
% 
% FILES NEEDED:
%  DUT_Runsheet.xlsx

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_Regimes

clear all
close all
clc


% READ RAFTING AND COLLISIONS DATA
dat = xlsread('DUT_Runsheet','E122:X177');

lambda_nd = dat(:,7);
H_target = dat(:,8);
H = dat(:,9);
sep = dat(:,12);
flag_col = dat(:,14);
flag_raf = dat(:,15);



%% 
% PLOT COLLISIONS/RAFTING/NIL
%  WAVELENGTH* VS WAVE HEIGHT

figure(1)
hold on

fill([0.8 0.8 1.8 1.8],[85 30 30 85],'r')
fill([1.8 1.8 3.3 3.3],[85 40 65 85],'g')
fill([3.3 3.3 5 5.5],[85 65 65 85],'m')
fill([1.8 5.2 5.5 7.5 7.5 5.5 5 3.3 1.8],[30 30 45 50 85 85 65 65 40],'c')
fill([5.2 7.5 7.5 5.5],[30 30 50 45],'y')


for j = 1:length(lambda_nd)
 if sep(j) == 30
  if H_target(j) == 40
   if flag_col(j) == 1 && flag_raf(j) == 0
    h1 = plot(lambda_nd(j),40,'bs','MarkerSize',10);
   elseif flag_raf(j) == 1
    h2 = plot(lambda_nd(j),40,'b^','MarkerSize',10);
   else
    h3 = plot(lambda_nd(j),40,'bx','MarkerSize',15);
   end
  elseif H_target(j) == 60
   if flag_col(j) == 1 && flag_raf(j) == 0
    plot(lambda_nd(j),60,'bs','MarkerSize',10)
   elseif flag_raf(j) == 1
    plot(lambda_nd(j),60,'b^','MarkerSize',10)
   else
    plot(lambda_nd(j),60,'bx','MarkerSize',15)
   end
  elseif H_target(j) == 80
   if flag_col(j) == 1 && flag_raf(j) == 0
    plot(lambda_nd(j),80,'bs','MarkerSize',10)
   elseif flag_raf(j) == 1
    plot(lambda_nd(j),80,'b^','MarkerSize',10)
   else
    plot(lambda_nd(j),80,'bx','MarkerSize',15)
   end
  end
 end
end

xlim([0 8])
ylim([20 100])
xlabel('Wavelength / Floe Length')
ylabel('Wave Height [mm]')
title('Separation = 30mm')

legend([h1 h2 h3],{'Collisions','Rafting','No Contact'})


figure(2)
hold on

fill([0.8 0.8 1.3 1.3 1.7 1.7],[85 30 30 50 60 85],'r')
fill([1.3 2 2 3.2 3.2 1.7 1.7 1.3],[30 30 45 50 85 85 60 50],'g')
fill([3.2 3.2 5 5.5],[85 70 70 85],'m')
fill([3.2 3.2 4 5],[70 50 50 70],'c')
fill([2 2 7.5 7.5 5.5 5 4 3.2],[45 30 30 85 85 70 50 50],'y')

for j = 1:length(lambda_nd)
 if sep(j) == 60
  if H_target(j) == 40
   if flag_col(j) == 1 && flag_raf(j) == 0
    h1 = plot(lambda_nd(j),40,'bs','MarkerSize',10);
   elseif flag_raf(j) == 1
    h2 = plot(lambda_nd(j),40,'b^','MarkerSize',10);
   else
    h3 = plot(lambda_nd(j),40,'bx','MarkerSize',15);
   end
  elseif H_target(j) == 60
   if flag_col(j) == 1 && flag_raf(j) == 0
    plot(lambda_nd(j),60,'bs','MarkerSize',10)
   elseif flag_raf(j) == 1
    plot(lambda_nd(j),60,'b^','MarkerSize',10)
   else
    plot(lambda_nd(j),60,'bx','MarkerSize',15)
   end
  elseif H_target(j) == 80
   if flag_col(j) == 1 && flag_raf(j) == 0
    plot(lambda_nd(j),80,'bs','MarkerSize',10)
   elseif flag_raf(j) == 1
    plot(lambda_nd(j),80,'b^','MarkerSize',10)
   else
    plot(lambda_nd(j),80,'bx','MarkerSize',15)
   end
  end
 end
end

xlim([0 8])
ylim([20 100])
xlabel('Wavelength / Floe Length')
ylabel('Wave Height [mm]')
title('Separation = 60mm')
legend([h1 h2 h3],{'Collisions','Rafting','No Contact'})






