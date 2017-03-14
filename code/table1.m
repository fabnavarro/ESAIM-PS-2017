% table1 - produces Table 1. of the paper
% Slope heuristics and V-Fold model selection 
% in heteroscedastic regression using strongly localized bases
%
%   Copyright (c) 2017 Fabien Navarro and Adrien Saumard

clear all
close all

rand('state',0);
randn('state',0);

nlist  = [256 1024 4096];
J      = log2(nlist);
MC     = 1000;
Pn     = @(x0,x)(mean((x-x0).^2,2));

noise_level = 'Low';
signal_name = {'Wave1','Wave2','HeaviSine1','HeaviSine2',...
               'Doppler1','Doppler2','Spikes1','Spikes2'};

j0 = 0;
V  = 2;
filter = 'Daubechies';
nVanishMmt = 8;
qmf = MakeONFilter(filter,nVanishMmt);

for in=1:length(signal_name)
  for im =1:length(nlist)
    for it=1:MC 
      [X,Y,s_star] = gendatafun(nlist(im),signal_name{in},noise_level);
      wc = FWT_PO(Y,j0,qmf);
      D  = 2.^(1:J(im)-1);
      hat_s_m  = zeros(length(D),nlist(im));
      wc_hat_s_m = zeros(1,nlist(im));

      n1 = nlist(im)/2;
      wc_odd_lin  = zeros(1,n1);
      wc_even_lin = zeros(1,n1);
      Y_odd = Y(1:2:end);
      Y_even = Y(2:2:end);
      
      wc_odd  = FWT_PO(Y_odd,j0,qmf);
      wc_even = FWT_PO(Y_even,j0,qmf);
      hat_s_odd  = zeros(length(D),n1);
      hat_s_even = zeros(length(D),n1);
      
      for ii =1:length(D)
        wc_hat_s_m(1:2^ii)=wc(1:2^ii);
        hat_s_m(ii,:) = IWT_PO(wc_hat_s_m,j0,qmf);

        wc_odd_lin(1:2^(ii))  = wc_odd(1:2^(ii));
        wc_even_lin(1:2^(ii)) = wc_even(1:2^(ii));
        hat_s_odd(ii,:)  = IWT_PO(wc_odd_lin,j0,qmf);
        hat_s_even(ii,:) = IWT_PO(wc_even_lin,j0,qmf);
      end
      err = Pn(hat_s_m,repmat(s_star,length(D),1));
      hat_sigma = sqrt(sum((Y-hat_s_m(end,:)).^2)/(nlist(im)-nlist(im)/2));      
      pen_Cp = 2*hat_sigma.^2*D'/nlist(im);
      mat_Y = repmat(Y,length(D),1);
      Crit_Cp  = Pn(hat_s_m,mat_Y) + pen_Cp;
      
      [~,m_star] = min(err);
      [~,m_cp] = min(Crit_Cp);
      hat_s_m_star = hat_s_m(m_star,:);
      hat_s_m_Cp   = hat_s_m(m_cp,:);
    
      hat_m = sh(D,nlist(im),Y,hat_s_m);
      hat_s_m_SH = hat_s_m(nextpow2(hat_m),:);
      
      % 2-CV Nason
      bar_s_odd = 0.5.*(hat_s_odd(:,1:n1-1) + hat_s_odd(:,2:n1));
      bar_s_odd(:,n1) = hat_s_odd(:,1);
      bar_s_even = 0.5.*(hat_s_even(:,1:n1-1) + hat_s_even(:,2:n1));
      bar_s_even(:,n1) = hat_s_even(:,1);
      
      mat_Y_odd = repmat(Y_odd,length(D),1);
      mat_Y_even = repmat(Y_even,length(D),1);
      Crit_CV = sum( (mat_Y_odd-bar_s_even).^2+...
                     (mat_Y_even-bar_s_odd).^2,2 );
      [~,hat_m_2FCV] = min(Crit_CV);
      hat_s_m_2FCV = hat_s_m(hat_m_2FCV,:);
      
      % 2-CVpen
      Crit_pen2F = 1/nlist(im)*sum((hat_s_m-mat_Y).^2,2)+...
      1/(2*nlist(im))*(-sum((mat_Y_even- hat_s_even).^2,2)...
                       +sum((mat_Y_odd - bar_s_even).^2,2)...
                       -sum((mat_Y_odd - hat_s_odd).^2,2)...
                       +sum((mat_Y_even- bar_s_odd).^2,2) );              
      [~,hat_m_pen2F] = min(Crit_pen2F);
      hat_s_m_pen2F = hat_s_m(hat_m_pen2F,:);
      
      MMSE(in,im,it)      = Pn(hat_s_m_star,s_star);
      MMSEcp(in,im,it)    = Pn(hat_s_m_Cp,s_star);
      MMSEsh(in,im,it)    = Pn(hat_s_m_SH,s_star);     
      MMSEcv(in,im,it)    = Pn(hat_s_m_2FCV,s_star);
      MMSEcvpen(in,im,it) = Pn(hat_s_m_pen2F,s_star);
    end
  end
end

Cor_sh = mean(MMSEsh./MMSE,3);
Cor_cp = mean(MMSEcp./MMSE,3);
Cor_cv = mean(MMSEcv./MMSE,3);
Cor_cvpen = mean(MMSEcvpen./MMSE,3);

Cor_sh_std = std(MMSEsh./MMSE,0,3)./sqrt(MC);
Cor_cp_std = std(MMSEcp./MMSE,0,3)./sqrt(MC);
Cor_cv_std = std(MMSEcv./MMSE,0,3)./sqrt(MC);
Cor_cvpen_std = std(MMSEcvpen./MMSE,0,3)./sqrt(MC);

Cor_table = [Cor_sh,Cor_cp,Cor_cv,Cor_cvpen];
Cor_table_sd = [Cor_sh_std,Cor_cp_std,Cor_cv_std,Cor_cvpen_std];

table(Cor_sh,Cor_cp,Cor_cv,Cor_cvpen,...
      'RowNames',signal_name)
