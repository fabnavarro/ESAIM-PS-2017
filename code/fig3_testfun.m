% fig3TestFun - produces Figure 3. of the paper
% Slope heuristics and V-Fold model selection 
% in heteroscedastic regression using strongly localized bases
%
%   Copyright (c) 2017 Fabien Navarro and Adrien Saumard

clear all
close all

rand('state',0);
randn('state',0);

n = 4096;
J = log2(n);
j0 = 0;
filter = 'Daubechies';
qmf = MakeONFilter(filter,8);
signal_names = {'Wave1','HeaviSine1','Doppler1','Spikes1'};       
              
for in=1:length(signal_names)
  figure(in);
  set(gcf,'Name',signal_names{in},'NumberTitle','off');
  subplot(121);
  [X,~,s_star] = gendatafun(n,signal_names{in},'Low');
  plot(X,s_star,'k','LineWidth',2)
  xlim([0 1]);ylim([0 1])
  xlabel('x');ylabel('s_*');
  subplot(122)
  wctheo = FWT_PO(s_star,j0,qmf);
  PlotWaveCoeff(wctheo,j0,1);
  axis([0 1 (-J) (-j0+1)])
end
tilefigs