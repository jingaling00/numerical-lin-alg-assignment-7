
% CHOLESKI'S ALGORITHM 6.6
% To factor the positive definite n by n matrix A into LL**T,
% where L is lower triangular.
A = [16 12 12 -4 -12 16 -4 0 0 -8;
    12 10 10 -4 -8 14 -3 2 -3 -8;
    12 10 14 4 0 18 -3 4 -1 -10;
    -4 -4 4 34 18 -6 5 10 15 12;
    -12 -8 0 18 30 2 9 14 3 2;
    16 14 18 -6 2 36 6 14 -12 -14;
    -4 -3 -3 5 9 6 45 15 16 19;
    0 2 4 10 14 14 15 46 5 4;
    0 -3 -1 15 3 -12 16 5 65 15;
    -8 -8 -10 12 2 -14 19 4 15 30];
b = [-20; -17; -17; 15; -1; -44; -32; -55; -43; 7];
 
 disp('tril(A)='), disp(single(tril(A)))
 pause
 N=length(A);
% STEP 1
 if A(1,1)<=0
    disp('A not positive definite')
    return
 end    
 A(1,1) = sqrt(A(1,1));
 single(tril(A))
 pause
% STEP 2
 for J = 2 : N 
 A(J,1) = A(J,1)/A(1,1);
 single(tril(A))
 pause
 end;
% STEP 3
 NN = N-1;
 for I = 2 : NN 
% STEP 4
 KK = I-1;
 S = A(I,1)*conj(A(I,1));
 for K = 2 : KK 
 S = S+A(I,K)*conj(A(I,K));
 end;
 if A(I,I)<=S
    disp('A not positive definite')
    return
 end    
 A(I,I) = sqrt(A(I,I)-S);
 single(tril(A))
 pause
% STEP 5
 JJ = I+1;
 for J = JJ : N 
 KK = I-1;
 S=A(J,1)*conj(A(I,1));
 for K = 2 : KK 
 S = S + A(J,K)*conj(A(I,K));
 end;
 A(J,I) = (A(J,I)-S)/A(I,I);
 single(tril(A))
 pause
 end;
 end;


% STEP 6
 S = A(N,1)*conj(A(N,1));
 for K = 2 : NN 
 S = S+A(N,K)*conj(A(N,K));
 end;
 if A(N,N)<=S
    disp('A not positive definite')
    return
 end   
 A(N,N) = sqrt(A(N,N)-S);
 single(tril(A)) 
 pause
 
 L=tril(A); 
 disp('L='), disp(single(L))

 % ----- MODIFICATIONS -----
% forward substitution for Ly = b
[n,~] = size(A);
y = [];
y(1) = b(1) / (L(1,1)); 
for i=2:n
    s = 0;
    for j=1:(i-1)
        s = s + (L(i,j) * y(j));
    end
    y(i) = (b(i) - s) / (L(i,i));
end

% back substitution for L^Tx = y
Lt = L.';
x = [];
x(n) = y(n) / (L(n,n)); 

for i=2:n
    k = n-i+1; % reverse mapping
    s = 0;
    for j=(k+1):n
        s = s + (Lt(k,j) * x(j));
    end
    fprintf('k %i s %i', k ,s );

    x(k) = (y(k) - s) / (L(k,k)); 
end
% END MODIFICATIONS


 
 
 
