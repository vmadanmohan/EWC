function TE = TransferEntFull(data,N,delay,teCalc)
    TE = zeros(N,N);
    for i = 1:N
        for j = 1:N
            if j~=i && ~isnan(delay(i,j))
                teCalc.initialise();
                D = num2str(delay(i,j));
                teCalc.setProperty('DELAY_PROP_NAME',D)
                teCalc.setObservations(data(:,i),data(:,j));
                I = teCalc.computeAverageLocalOfObservations();
                if I<0 || teCalc.computeSignificance().pValue>=0.01
                    TE(i,j)=0;
                else
                    TE(i,j)=I;
                end
            else
                TE(i,j)=0;
            end
        end
    end
