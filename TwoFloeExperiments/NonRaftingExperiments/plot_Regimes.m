% function plot_Regimes
%
% LJ YIEW
% Created on  May 2016
% Last edited Dec 2016
%
% Plots the collision regimes (according to number of collision events) as
% a function of wavelength and wave height
%
% FILES NEEDED:
%  wavefield.m
%  plotxx.m
%  Data_ExperimentVsModel.xlsx

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_Regimes

dat = xlsread('Data_ExperimentVsModel.xlsx','B6:H67');


freq = dat(:,2);
wh = dat(:,3);
wh_target = dat(:,4);
col = dat(:,5);

for j = 1:length(freq)
 if isnan(freq(j))
  freq(j) = 0;
 end
 [field] = wavefield('f',freq(j),0.5);
 lambda(j) = cell2mat(field(4,2)); % corresponding wavelength
end
lambda_nd = lambda./0.4;


for j = 1:length(col)
 if col(j) == 0
  col_flag(j) = 0;
  col_few_flag(j) = NaN;
 elseif isnan(col(j))
  col_flag(j) = NaN;
  col_few_flag(j) = NaN;
 else
  col_flag(j) = 1;
  if col(j) <= 3
   col_few_flag(j) = 1;
  else
   col_few_flag(j) = 0;
  end
 end
end


%%
% ARRANGE ACCORDING TO WAVELENGTH
figure
hold on
for j = 1:length(freq)
 if col_flag(j) == 1
  if col_few_flag(j) == 1
   h1 = plot(lambda_nd(j),wh(j),'bs','MarkerSize',10);
  elseif col_few_flag(j) == 0
   h2 = plot(lambda_nd(j),wh(j),'bo','MarkerSize',10);
  end
 elseif col_flag(j) == 0
  h3 = plot(lambda_nd(j),wh(j),'bx','MarkerSize',15);
 end
end
xlabel('Wavelength / Floe Length')
ylabel('Wave Height [mm]')
ylim([0 100])
xlim([0 12])
legend([h1 h2 h3],{'Collisions (\leq3)','Collisions (>3)','No Collisions'})

%%
% ARRANGE ACCORDING TO WAVE FREQUENCY
figure
hold on
for j = 1:length(freq)
 if col_flag(j) == 1
  if col_few_flag(j) == 1
   h1 = plot(freq(j),wh(j),'ks','MarkerSize',10);
  elseif col_few_flag(j) == 0
   h2 = plot(freq(j),wh(j),'ko','MarkerSize',10);
  end
 elseif col_flag(j) == 0
  h3 = plot(freq(j),wh(j),'kx','MarkerSize',15);
 end
end
xlabel('Wavelength / Floe Length')
ylabel('Wave Height [mm]')
ylim([0 100])
xlim([0.4 1.6])
legend([h1 h2 h3],{'Collisions (\leq3)','Collisions (>3)','No Collisions'})

ax1 = gca;
set(ax1,'XTick',[0.4 0.6 0.8 1 1.2 1.4 1.6])
set(ax1,'XTickLabel',{'0.4','0.6','0.8','1','1.2','1.4','1.6'})


% SECOND X SCALE FOR NONDIMENSIONAL WAVELENGTH
xlabels{1} = 'Wave Frequency [Hz]';
xlabels{2} = 'Wavelength/Floe Diameter';
ylabels{1} = 'Wave Height [mm]';
ylabels{2} = '';
hold on
plotxx(0,0,0,0,xlabels,ylabels);
xlim([0.4 1.6])
ylim([0 0.2])
ax = gca;
set(ax,'XTick',[0.4 0.6 0.8 1 1.2 1.4 1.6])
set(ax,'YTick',[0 0.05 0.1 0.15 0.2])
set(ax,'XTickLabel',{'16.2','9.5','5.9','3.9','2.7','2','1.5'})
set(ax,'YTickLabel',{'','','','',''})

return

