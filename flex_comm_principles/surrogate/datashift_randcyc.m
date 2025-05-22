% Circularly shifts data relative to a row 'k', by 'p' steps
function randcyc=datashift_randcyc(a,k,p)
    b = a(:,k);
    a(:,k)=[];
    a = circshift(a,p,1);
    randcyc = [a(:,1:k-1) b a(:,k:end)];
