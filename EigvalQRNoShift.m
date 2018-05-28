function [eigenvalues, i, s] = EigvalQRNoShift(D,tol, imax)
%tol - tolerancja (górna granica wartosci) elementów zerowanych
%imax - maksymalna liczba iteracji
n=size(D, 1);
i=1;
while i <= imax && max(max(D-diag(diag(D)))) > tol
    [Q1, R1]=qr(D);
    D=R1*Q1; %macierz przekszta?cona
    i=i+1;
end
if i >=imax
        s=0;
    else
        s=1;
end
if i > imax
    error('imax exceeded program terminated');
end
eigenvalues=diag(D);
end

