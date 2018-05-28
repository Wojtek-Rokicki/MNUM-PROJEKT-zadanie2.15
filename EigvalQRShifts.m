function [eigenvalues, iteracje, s] = EigvalQRShifts(A, tol, imax)
%tol - tolerancja (gorna granica wartosci) elementow zerowych
%imax - max liczba iteracji dla liczenia jednej wartosci wlasnej
n=size(A,1);
eigenvalues=diag(zeros(n));
INITIALsubmatrix=A; %macierz poczatkowa (oryginalna)
iteracje=0;
for k=n:-1:2
    DK=INITIALsubmatrix; %macierz startowa dla jednej wartosci wlasnej
    i=0;
    while i <= imax && max(abs(DK(k,1:k-1))) > tol
        DD=DK(k-1:k,k-1:k); %2x2 podmacierz prawego dolnego rogu
        [ev1,ev2]=quadpolynroots(1,-(DD(1,1)+DD(2,2)),DD(2,2)*DD(1,1)-DD(2,1)*DD(1,2));
        if abs(ev1-DD(2,2)) < abs(ev2-DD(2,2))
            shift=ev1; %najblizsza DK(k,k) wartosc wlasna podmacierzy DD
        else
            shift=ev2; %najblizsza DK(k,k) wartosc wlasna podmacierzy DD
        end
        DK=DK-eye(k)*shift; %macierz przesunieta
        [Q1,R1]=qr(DK); %faktoryzacja QR
        DK=R1*Q1+eye(k)*shift; %macierz przeksztalcona
        i=i+1;
        iteracje = iteracje + 1;
    end
    if i >=imax
        s=0;
    else
        s=1;
    end
    if i > imax
        error('imax exceeded program terminated');
    end
    eigenvalues(k)=DK(k,k);
    if k >2
        INITIALsubmatrix=DK(1:k-1,1:k-1); %deflacja macierzy
    else
        eigenvalues(1)=DK(1,1); %ostatnia wartosc wlasna
    end
end

end

