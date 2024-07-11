function [cMI,cMI_std] = CondMutInfo(data,win,N,delay,cmiCalc)
    thresh=3;           % Z-score threshold
    cMI = zeros(N,N);
    cMI_std = zeros(N,N);
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
    for i=1:N          % Loop over source regions
        reldelays = delay(i,:);  % Delays relative to the source region
        delcorr = zeros(max(reldelays,[],'omitnan')+ size(data,1),N);      % Delay corrected matrix
        signum = sum(sigevents(:,i),1);  % Number of significant events
        if signum~=0
            sigind = find(sigevents(:,i))+max(reldelays,[],'omitnan');  % Indices of significant events
            totevents=signum;
            cMI_event = zeros(totevents,N);
            for k=1:N                                   % Shifting time series based on delays
                if k~=i && ~isnan(reldelays(k))  
                    delcorr(end-reldelays(k)-size(data,1)+1:end-reldelays(k),k)=data(:,k);
                else
                    delcorr(end-size(data,1)+1:end,k) = data(:,k);
                end
            end
            while signum            % While loop starting from the last significant event and moving backwards
                reorg = zeros(win,3);
                reorg(:,1) = delcorr(sigind(signum):sigind(signum)+win-1,i);
                for j=1:N   % Loop over all possible targets
                    if j~=i && ~isnan(reldelays(j))
                        reorg(:,2) = delcorr(sigind(signum):sigind(signum)+win-1,j);
                        reorg(:,3) = delcorr(sigind(signum)-win:sigind(signum)-1,j);% Conditioning on past of target region (in this case, j)
                        cmiCalc.initialise(1,1,1);
                        cmiCalc.setObservations(reorg(:,1),reorg(:,2),reorg(:,3));
                        I = cmiCalc.computeAverageLocalOfObservations();  % (Conditional) Mutual Information calculation
                        if I<0 || cmiCalc.computeSignificance().pValue>=0.01/totevents  % Applying a Bonferroni correction to account for multiple events
                            cMI_event(signum,j)=0;
                        else
                            cMI_event(signum,j)=I;
                        end
                    else
                        cMI_event(signum,j)=0;
                    end
                end
                signum=signum-1;
            end
            cMI(i,:)=mean(cMI_event,1);  % Mean cMI over all significant events in source region
            cMI_std(i,:)=std(cMI_event,0,1);
        end
    end
