% Function to identify significant communication events using a z-score threshold. The locations of the events are used to estimate EWC and neural oscillatory measures in both "observed" and "surrogate" datasets
function sigevents = eventiden(data,delay,win,thresh)
    sigevents = abs(zscore(data))>=thresh;
    sigevents(1:win,:)=0;                       % Remove sig. events 1 window at start of data (to make room for the self-conditioning window)
    sigevents(end-win-max(max(delay,[],'omitnan')):end,:)=0; % Remove sig. events 1 window + max(delay) at the end of the data, to make room for windows and to avoid computations corrupted by the zero padding
    for i=1:size(sigevents,2)                   % Remove additional significant events that fall inside a window
        for j=1:size(sigevents,1)
            if sigevents(j,i)==1
                sigevents(j+1:j+win+1,i)=0;
                j=j+win+1;
            end
        end
    end
