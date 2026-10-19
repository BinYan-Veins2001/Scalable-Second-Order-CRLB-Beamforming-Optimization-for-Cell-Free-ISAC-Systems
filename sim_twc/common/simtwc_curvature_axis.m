function y=simtwc_curvature_axis(x,transition,ax)
%SIMTWC_CURVATURE_AXIS Signed-log display with raw-value tick labels.

y=sign(x).*log10(1+abs(x)/transition);
if nargin<3 || isempty(ax); return; end

limit=ylim(ax);
positive=transition*[10 100 1e4 1e6 1e8 1e10 1e12];
values=[-fliplr(positive) 0 positive];
positions=sign(values).*log10(1+abs(values)/transition);
use=positions>=limit(1) & positions<=limit(2);
values=values(use); positions=positions(use);
labels=cell(size(values));
for k=1:numel(values)
    if values(k)==0
        labels{k}='0';
    else
        exponent=round(log10(abs(values(k))));
        if exponent==0
            labels{k}=sprintf('%s1',sign_text(values(k)));
        else
            labels{k}=sprintf('%s10^{%d}',sign_text(values(k)),exponent);
        end
    end
end
set(ax,'YTick',positions,'YTickLabel',labels,'TickLabelInterpreter','tex');
end

function s=sign_text(x)
if x<0; s='-'; else; s=''; end
end
