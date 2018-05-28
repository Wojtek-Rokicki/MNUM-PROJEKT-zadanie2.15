function [ w ] = pval(a , x)
%   pval oblicza wartosci wielomianu o wspolczynnikach a(a(1) odpowiada x^0)
%   w - wartosci wielomianu w danych punktach x
    ilprobek = size(x,1);
    stwiel = size(a,1);
    w = zeros(ilprobek,1);
    for i = 1: ilprobek
        for j = 1:stwiel
            w(i) = w(i) + a(j) * x(i)^(j-1);
        end
    end
end
