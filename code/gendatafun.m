function [X, Y, sig] = gendatafun(n,signal_name,noise_level)
% Generates a data sample of size n corresponding to experiment.
%
%%%
% This function is based on MakeSinalNewb written by
% Anestis Antoniadis and Jeremie Bigot
% which is based on a code provided to them by Buckheit, Chen, Donoho,
% Johnstone & Scargle.
%%%
%
% 	Input
%		n	        Desired sampling length.
%		signal_name	String: 'Wave1|2','HeaviSine1|2','Doppler1|2','Spikes1|2'.
%		noise_level Level of noise.
%
%	Output
%		X   sorted X_i.
%		Y   observed data Y_i.
%       sig s_star.
%
%
%   Copyright (c) 2017 Fabien Navarro and Adrien Saumard


X = sort(rand(1,n));
epsilon = randn(1,n);
if strcmp(signal_name,'Step1') || strcmp(signal_name,'Step2'),
  sig = 0.2 + 0.6*(X > 1/3 & X <= 0.75); 
elseif strcmp(signal_name,'Wave1') || strcmp(signal_name,'Wave2'),
  sig = 0.5 + (0.2.*cos(4*pi*X)) + (0.1.*cos(24*pi*X));
elseif strcmp(signal_name,'HeaviSine1')|| strcmp(signal_name,'HeaviSine2'),
  sig = 4.*sin(4*pi.*X) - sign(X - .3) - sign(.72 - X) + 5;
  sig = (0.6/9)*sig + 0.2;
elseif strcmp(signal_name,'Doppler1')|| strcmp(signal_name,'Doppler2'),
  sig = sqrt(X.*(1-X)).*sin((2*pi*1.05) ./(X+.05)) + 0.5;
  sig = 0.6*sig + 0.2;
elseif strcmp(signal_name,'Spikes1')|| strcmp(signal_name,'Spikes2'),
  sig = 15.6676 .* (exp(-500.*(X-0.23).^2) + ...
		2.*exp(-2000.*(X-0.33).^2) +  4.*exp(-8000.*(X-0.47).^2) + ...
		3.*exp(-16000.*(X-0.69).^2) +exp(-32000.*(X-0.83).^2) );
		sig=(0.6/range(sig)).*sig+0.2;		
else
  disp('Unknown signal names');
  disp('Allowable Names are:')
  disp('Wave1|2'),
  disp('HeaviSine1|2'),
  disp('Doppler1|2'),
  disp('Spikes1|2')
end
 
if strcmp(noise_level,'Low')
  kap = 0.01;
elseif strcmp(noise_level,'High')
  kap = 0.05;
else
  disp('Unknown noise level')
  disp('Allowable Names are:')
  disp('Low'),
  disp('High')
end
    
switch signal_name
  case {'Step1','Wave1','Blip1','Blocks1','Bumps1',...
		'HeaviSine1','Doppler1','Angles1',...
	    'Parabolas1','TShSine1','Spikes1','Corner1'}
    sigma_v = kap;
  case {'Step2','Wave2','Blip2','Blocks2','Bumps2',...
		'HeaviSine2','Doppler2','Angles2',...
	    'Parabolas2','TShSine2','Spikes2','Corner2'}
    sigma_v = X*kap*2;
  otherwise
     sigma_v=NaN*ones(1,n);
end

Y = sig+sigma_v.*epsilon;
end
