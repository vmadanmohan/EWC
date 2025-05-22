function [P,P_std] = PearsonEWCNetwork_cyclicsurr(orig_data,SC,win,N,delay,sigevents,r)
    P = zeros(N,N);
    P_std = zeros(N,N);
    for i=1:N          % Loop over source regions
        data=datashift_randcyc(orig_data,i,r(i));
	reldelays = delay(i,:);  % Delays relative to the source region
        signum = sum(sigevents(:,i),1);  % Number of significant events
        if signum~=0
            sigind = find(sigevents(:,i));  % Indices of significant events
            totevents=signum;
            P_event = zeros(totevents,N);
            while signum            % While loop starting from the last significant event and moving backwards
                for j=find(SC(i,:))   % Loop over regions that are anatomically connected to the source
                    [I,pval] = partialcorr(data(sigind(signum):sigind(signum)+win-1,i),data(sigind(signum)+reldelays(j):sigind(signum)+reldelays(j)+win-1,j),data(sigind(signum)+reldelays(j)-win:sigind(signum)+reldelays(j)-1,j));  % Partial correlation calculation
                    if pval>=0.01/totevents  % Applying a Bonferroni correction to account for multiple events
                        P_event(signum,j)=0;
                    else
                        P_event(signum,j)=I;
                    end
                end
                signum=signum-1;
            end
            P(i,:)=mean(P_event,1);  % Mean PC over all significant events in source region
            P_std(i,:)=std(P_event,0,1);
        end
    end
