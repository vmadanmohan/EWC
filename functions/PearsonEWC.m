function [P,P_std] = PearsonEWC(data,win,N,delay)
    thresh=3;           % Z-score threshold
    P = zeros(N,N);
    P_std = zeros(N,N);
    sigevents = abs(zscore(data))>=thresh;        % Event identification = z-scoring and thresholding
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
            P_event = zeros(totevents,N);
            for k=1:N                                   % Shifting time series based on delays
                if k~=i && ~isnan(reldelays(k))  
                    delcorr(end-reldelays(k)-size(data,1)+1:end-reldelays(k),k)=data(:,k);
                else
                    delcorr(end-size(data,1)+1:end,k) = data(:,k);
                end
            end
            while signum            % While loop starting from the last significant event and moving backwards
                for j=1:N   % Loop over all possible targets
                    if j~=i && ~isnan(reldelays(j))
                        [I,pval] = partialcorr(delcorr(sigind(signum):sigind(signum)+win-1,i),delcorr(sigind(signum):sigind(signum)+win-1,j),delcorr(sigind(signum)-win:sigind(signum)-1,j));  % Partial correlation calculation
                        if pval>=0.01/totevents  % Applying a Bonferroni correction to account for multiple events
                            P_event(signum,j)=0;
                        else
                            P_event(signum,j)=I;
                        end
                    else
                        P_event(signum,j)=0;
                    end
                end
                signum=signum-1;
            end
            P(i,:)=mean(P_event,1);  % Mean PC over all significant events in source region
            P_std(i,:)=std(P_event,0,1);
        end
    end
