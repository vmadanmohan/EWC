%% Left Hemisphere
c=[];
for sub=1:30
    T=TE_mean(1:50,1:50,sub);
    P=sqrt(P_mean(1:50,1:50,sub).^2);
    c(sub)=corr(P(T>0),T(T>0),"type","Pearson");
end

figure;histogram(c,7,"Normalization","pdf");xline(median(c));hold on;ksdensity(c,'Bandwidth',0.02)
xlabel("|Pearson-R|")
ylabel("Probability Density")
legend("","Median correlation")
fontsize(gcf,18,"points")

P=sqrt(P_mean(1:50,1:50,16).^2);    % Taking SUB 16 as a representative since correlation value closest to median
T=TE_mean(1:50,1:50,16);
figure;scatter(P(T>0),T(T>0));lsline
y=fitlm(P(T>0),T(T>0));
xlabel("Partial corr. (EWC)")
ylabel("TE (Full)")
legend('',sprintf('%.4fx+%.4f',y.Coefficients.Estimate(2),y.Coefficients.Estimate(1)))
annotation('textbox',"String",sprintf("R^2=%.4f",y.Rsquared.Ordinary))
fontsize(gcf,18,"points")
[R,rho,RLO,RUP]=corrcoef(P(T>0),T(T>0))

figure;boxchart([(runtimeP./numepochs)/2,(runtimeTE./numepochs)/2]);
ylabel("Inference time per epoch (s)")
xticklabels(["Partial corr. (EWC)","TE (Full)"])
fontsize(gcf,18,"points")

%% Right Hemisphere
c=[];
for sub=1:30
    T=TE_mean(51:end,51:end,sub);
    P=sqrt(P_mean(51:end,51:end,sub).^2);
    c(sub)=corr(P(T>0),T(T>0),"type","Pearson");
end

figure;histogram(c,7,"Normalization","pdf");xline(median(c));hold on;ksdensity(c,'Bandwidth',0.01)
xlabel("|Pearson-R|")
ylabel("Probability Density")
legend("","Median correlation")
fontsize(gcf,18,"points")

P=sqrt(P_mean(51:end,51:end,11).^2);    % Taking SUB 11 as a representative since correlation value closest to median
T=TE_mean(51:end,51:end,11);
figure;scatter(P(T>0),T(T>0));lsline
y=fitlm(P(T>0),T(T>0));
xlabel("Partial corr. (EWC)")
ylabel("TE (Full)")
legend('',sprintf('%.4fx+%.4f',y.Coefficients.Estimate(2),y.Coefficients.Estimate(1)))
annotation('textbox',"String",sprintf("R^2=%.4f",y.Rsquared.Ordinary))
fontsize(gcf,18,"points")
[R,rho,RLO,RUP]=corrcoef(P(T>0),T(T>0))

figure;boxchart([(runtimeP./numepochs)/2,(runtimeTE./numepochs)/2]);
ylabel("Inference time per epoch (s)")
xticklabels(["Partial corr. (EWC)","TE (Full)"])
fontsize(gcf,18,"points")
