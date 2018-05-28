function [] = Zadanie1(rozm, ilosc)
imax = 1000;
tol = 0.00001;
isumsa=0; %suma iteracji z algorytmow z przesunieciem niesymetrycznej macierzy
isumns=0; %suma iteracji z algorytmow bez przesuniecia
isums=0; %suma iteracji z algorytmow z przesunieciem
failsa=0;
failns=0;
fails=0;
for i=1:ilosc
    A=rand(rozm);
    
    [e, iter, s]=EigvalQRShifts(A, tol, imax);
    isumsa=isumsa+iter;
    if s==0
        failsa=failsa+1;
    end
        
    disp(e);
    
    A=A+A';
    [e, iter, s]=EigvalQRNoShift(A, tol, imax);
    isumns=isumns+iter;
    if s==0
        failns=failns+1;
    end
        
    disp(e);
    [e, iter, s]=EigvalQRShifts(A, tol, imax);
    isums=isums+iter;
    if s==0
        fails=fails+1;
    end
    disp(e);
    e=eig(A);
    disp(e);
end

disp(isumsa/ilosc);
disp(failsa);
disp(isumns/ilosc);
disp(failns);
disp(isums/ilosc);
disp(fails);
    
end