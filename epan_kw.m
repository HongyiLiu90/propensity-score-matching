function [oc0 oc1] = epan_kw(h,Y,treat,p)

% Epanechnikov kernel weight routine
% calculates PS estimate of E[y|x] at each point 
% given in the vector w.
%
%
% Y: an n x 1 vector
% w: an  n x 1 vector
% oc1: an n x 1 vector 
% oc1: an n x 1 vector
% note this program doesn't check to make sure inputs fit these 
% requirements

n = length(Y);
[~,order]=sort(treat);
sp = p(order,:); % sort the data matrix by treated variable
sY = Y(order,:); % sort the data matrix by treated variable
streat = treat(order,:);
c = find(streat==0);
n0 = length(c);
i1 = n0 + 1;

w = zeros(n,1);
k =[];
dist =[];

% Matching each observation in group=1
for i=i1:n,
   
% construct standartized distance and compute kernel weights among group=0
   distance = abs( sp-sp(i,1) )/h;
   dist = distance(streat==0);
   k =(dist < 1).*(1-dist.^2)*3/4;
   
   if mean(k)== 0
       w(i,1) = 0;
   else
       w(streat == 0) = w(streat == 0) + k/(mean(k)*n0);
       w(i,1) = 1;
   end
  
     
end;
   % calculate an index based on distance to w
   wtilde = w(streat==0)/sum(w(streat==0));
   oc0 = sY(streat == 0).*wtilde*n0;
   % mean outcome among the non-treated matched to the treated and weighted
   % by the distribution of ps among treated
   oc1 = sY((streat ==1)&(w~=0));
   % mean outcome among the treated for whom at least one similar 
   % non-treated is found

   
   
   
   
   