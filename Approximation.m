function [ a, res ] = Approximation( x, y, n, meth )
%   n - stopien wielomianu
%   meth : 1 - uklad rownan normalnych; 2 - qr
    N = size(x,1);
    A = zeros(N,n);
    
    %Wypelniamy macierz A odpowiednimi potêgami x
    for i=1:N
        for j = 1:n
            A(i,j) = x(i,1)^(j-1);
        end
    end
    
    %uklad rownan normalnych
    if meth == 1
        ata = A'*A;
        aty = A'*y;
        a = ata\aty;
        res = norm(aty - ata*a,Inf);
    %uklad wynikajacy z rozkladu QR
    elseif meth == 2
        [Q,R] = QRDecomposition(A);
        a =R\Q'*y;
        res = norm(R*a - Q'*y,Inf);
    end
end

