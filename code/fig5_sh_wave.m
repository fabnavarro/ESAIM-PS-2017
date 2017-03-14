% fig5SH - produces Figure 5. of the paper
% Slope heuristics and V-Fold model selection 
% in heteroscedastic regression using strongly localized bases
%
% Copyright (c) 2017 Fabien Navarro and Adrien Saumard

clear all
close all

signal_name = {'Wave1','Wave2'};
noise_level = {'Low','High'};

Pn = @(x0,x)(mean((x-x0).^2,2));

n = 1024;
J = log2(n);
j0 = 0;
D = 2.^(1:J-1);
filter = 'Daubechies';
qmf = MakeONFilter(filter,8);

for im=1:length(noise_level)
  for in=1:length(signal_name)
    rand('seed',0)
    randn('seed',0)
    hat_s_m = zeros(length(D),n);
    wc_hat_s_m = zeros(1,n);
    [X,Y,s_star] = gendatafun(n,signal_name{in},noise_level{im});
    wc = FWT_PO(Y,j0,qmf);
    for ii =1:length(D)
      wc_hat_s_m(1:2^ii) = wc(1:2^ii);
      hat_s_m(ii,:) = IWT_PO(wc_hat_s_m,j0,qmf);
    end
    err = Pn(hat_s_m,repmat(s_star,length(D),1));
    [~,m_star] = min(err);
    hat_s_m_star = hat_s_m(m_star,:);
    [m_sh,critSH] = sh(D,n,Y,hat_s_m);
    m_sh = nextpow2(m_sh);
    hat_s_sh = hat_s_m(m_sh,:);

    figure;
    set(gcf,'Name',[signal_name{in} noise_level{im}],'NumberTitle','off');
    subplot(221);    
    plot(X,Y,'k','LineWidth',2);
    xlim([0 1]);ylim([0 1]);
    xlabel('$X$','Interpreter','latex');
    ylabel('$Y$','Interpreter','latex');
    subplot(222);
    plot(X,s_star,'k--','LineWidth',2);hold on
    plot(X,hat_s_sh,'k','LineWidth',2);
    legend('$s_*$','$\widehat{s}_{\widehat{m}_\mathrm{SH}}$',...
           'Location','NorthEast','Orientation','Vertical');
    legend boxoff;
    h0 = legend;set(h0, 'interpreter', 'latex','FontSize',15)
    xlim([0 1]);ylim([0 1]);
    subplot(223)
    rsh = rescale(critSH,min(err),max(err));
    loglog(D,err,'k-o','LineWidth',1.2);grid;hold on;axis tight
    h = loglog(D,rsh,'k-o','LineWidth',1.2);
    set(h, 'color', [0.5 0.5 0.5])
    loglog(D(m_sh),rsh(m_sh),'o','LineWidth',2,...
           'MarkerEdgeColor','k','MarkerFaceColor',[.75 .75 .75],...
           'MarkerSize',12);
    loglog(D(m_star),err(m_star),'k*','LineWidth',2,'MarkerSize',6)
    xlabel('$D_m$', 'interpreter', 'latex');
    legend('$\ell(s_{\ast},\widehat{s}_{m})$',...
           '$\mathrm{crit_{\mathrm{SH}}}(m)$',...
           ['$D_{\widehat{m}_\mathrm{SH}}=$' num2str(2^m_sh)],...
           ['$D_{m_{\ast}}=$' num2str(2^m_star)],...
           'Location','NorthEast','Orientation','Vertical');
    h1 = legend;set(h1, 'interpreter', 'latex','FontSize',12)
    subplot(224)
    wcn = wc_hat_s_m;
    wcn(D(m_sh)+1:end) = 0;
    PlotWaveCoeff(wcn,j0,1);
  end
end
tilefigs
