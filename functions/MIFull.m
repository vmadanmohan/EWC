function MI = MIFull(data,N,delay,miCalc)
    MI=zeros(N,N);
    for i = 1:N
        for j = 1:N
            if j~=i && ~isnan(delay(i,j))
                miCalc.initialise(1,1);
                miCalc.setObservations(data(:,i),data(:,j));
                I = miCalc.computeAverageLocalOfObservations();
                if I<0 || miCalc.computeSignificance().pValue>=0.01
                    MI(i,j)=0;
                else
                    MI(i,j)=I;
                end
            else
                MI(i,j)=0;
            end
        end
    end
