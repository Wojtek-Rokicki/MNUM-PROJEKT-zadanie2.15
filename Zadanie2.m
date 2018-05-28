function [ r, max ] = Zadanie2( n , meth)
    if meth == 1
        rodz = 'n';
        tytul  = ['Uklad rownan normalnych - stopien wielomainu = ' int2str(n-1)];
    else
        rodz = 'qr';
        tytul  = ['Rozklad QR - stopien wielomainu = ' int2str(n-1)];
    end
    gx = -5:0.1:5;
    gx = gx';
    [x,y] = data();
    [a,r] = Approximation(x,y,n,meth);
    z = pval(a,gx);
    max=0;
    i=1;
    j=1;
    while i~=11
        if abs(y(i)-z(j))>max
            max=abs(y(i)-z(j));
        end
        i=i+1;
        j=j+10;
    end
    scatter(x,y);
    hold on;
    plot(gx,z);
    title(tytul);
    hold off;

end