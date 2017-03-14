function [Modelselectmax,critSH] = sh(D,n,Y_v,flin)
% The slope heuristics method based on dimension jump approach.
%
%%%
% This function is based on CAPUSHE
% Authors:J.P. Baudry(Paris XI), C.Maugis(IMT) and B. Michel(LSTA)
% email: cathy.maugis-at-insa-toulouse.fr
% Web: http://www.math.univ-toulouse.fr/~maugis/CAPUSHE.html
% February 2010
%%%
%
% 	Input
%		D	 models dimension.
%		n	 sample size.
%		Y_v	 observed data.
%		flin models collection.
%
%	Output
%		Modelselectmax	selected model.
%		critSH          slope heuristics criterion.
%

min_cont   = sum((repmat(Y_v,length(D),1)-flin).^2,2)/n;
ModelContr = min_cont';
ModelShape = D/n;
ModelComplexity = D;
ModelNames = D;

KappaPath=[];Modelpath=[];

%%Determination of the jumps
index=1:length(ModelContr);
%Initialization
i=1;
KappaPath(i)=0;
Modelpath(i)=argmin(index,ModelContr+KappaPath(i)*ModelShape,ModelShape);

%Step i, i > 1:
while isfinite(KappaPath(i))
i=i+1;
Modelprev=Modelpath(i-1);
G=index(( ModelContr > ModelContr(Modelprev)) & (ModelShape < ModelShape(Modelprev)));
if length(G)>0
  KappaPath(i)=min((ModelContr(G) - ModelContr(Modelprev))./(ModelShape(Modelprev) - ModelShape(G)));
  Modelpath(i)=argmin(G,ModelContr(G)+KappaPath(i)*ModelShape(G),ModelShape(G));
else
KappaPath(i)=Inf;
end
end
Complexpath=ModelComplexity(Modelpath)';

%%Selection of the maximal jump
jumps=Complexpath(1:(end-1)) - Complexpath(2:end);
ijumps=find(jumps==max(jumps));
KappaSelectMax=KappaPath(max(ijumps)+1);
critSH = ModelContr + (2*KappaSelectMax*ModelShape);
[bb,b3]=min(critSH);
Modelselectmax=ModelNames(b3);
end


function x=argmin(l1,l2,l3)
% This function returns the value of l1 corresponding to the minimum of l2
% If l2 has several minimas, take the one for which l3 is minimal
% If there are still several possible values, take the 'min' of them
%
% (C) Sylvain Arlot, 2009
%
% INPUT:
% 3 vectors l1, l2, l3 of same size
%
% OUTPUT:
% x : real number = value of l1 corresponding to the minimum of l2
%
index=(1:numel(l1));
x1=index(l2==min(l2));% positions of the minimum of l2
if numel(x1)==1
    x=l1(x1);
else
    index1=(1:numel(x1));
    l31=l3(x1);
    l11=l1(x1);
    x2=index1(l31==min(l31));
    x=min(l11(x2));
end
end
