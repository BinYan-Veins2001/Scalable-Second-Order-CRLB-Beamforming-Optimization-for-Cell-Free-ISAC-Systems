function [mu,se,n] = simtwc_mean_se(x,dim)
%SIMTWC_MEAN_SE NaN-aware mean and standard error.

if nargin<2; dim=1; end
n=sum(isfinite(x),dim);
mu=mean(x,dim,'omitnan');
sd=std(x,0,dim,'omitnan');
se=sd./sqrt(max(n,1));
se(n==0)=NaN;
end
