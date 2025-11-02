% DIRECT FACTORIZATION ALGORITHM 6.4
% To factor the n by n matrix A as A=LU
% ISW=0 for Doolittle; ISW=1 for Crout 
 A = [2 4 0 0 4 2 2 2;
     4 10 2 0 10 8 6 6;
     0 2 4 2 6 4 4 6;
     0 0 2 3 6 1 2 3;
     4 10 6 6 19 7 10 15;
     2 8 4 1 7 10 8 12;
     2 6 4 2 10 8 10 16;
     2 6 6 3 15 12 16 27];
 disp('A='), disp(single(A))
 pause
 N=length(A);
 ISW = 0;
 TRUE = 1;
 FALSE = 0;
 OK=TRUE; 
 for I = 1 : N 
 XL(I) = 1;
 end;
% STEP 1
 if abs(A(1,1)) <= 1.0e-20 
 OK = FALSE;
 else
% the entries of L below the main diagonal will be placed 
% in the corresponding entries of A; the entries of U 
% above the main diagonal will be placed in the 
% corresponding entries of A; the main diagonal which 
% was not input will become the main diagonal of A; 
% the input main diagonal of L or U is placed in XL
 A(1,1) = A(1,1)/XL(1);
% STEP 2
 for J = 2 : N 
 if ISW == 0 
% first row of U
 A(1,J) = A(1,J)/XL(1);
% first column of L
 A(J,1) = A(J,1)/A(1,1);
 single(A)
 pause
 else
% first row of U
 A(1,J) = A(1,J)/A(1,1);
% first column of L
 A(J,1) = A(J,1)/XL(1);
 single(A)
 pause 
 end;
 end;
% STEP 3
 M = N-1;
 I = 2;
 while I <= M & OK == TRUE 
% STEP 4
 KK = I-1;
 S = A(I,1)*A(1,I);
 for K = 2 : KK 
 S = S+A(I,K)*A(K,I);
 end;
 A(I,I) = (A(I,I)-S)/XL(I);
 single(A)
 pause
 if abs(A(I,I)) <= 1.0e-20 
 OK = FALSE;
 else
% STEP 5
 JJ = I+1;
 for J = JJ : N 
 SS = A(I,1)*A(1,J);
 S = A(J,1)*A(1,I);
 for K = 2 : KK 
 SS = SS+A(I,K)*A(K,J);
 S = S+A(J,K)*A(K,I);
 end;
 if ISW == 0 
% Ith row of U
 A(I,J) = (A(I,J)-SS)/XL(I);
% Ith column of L
 A(J,I) = (A(J,I)-S)/A(I,I);
 single(A)
 pause
 else
% Ith row of U
 A(I,J) = (A(I,J)-SS)/A(I,I);
% Ith column of L
 A(J,I) = (A(J,I)-S)/XL(I);
 single(A)
 pause 
 end;
 end;
 end;
 I = I+1;
 end;
 if OK == TRUE 
% STEP 6
 S = A(N,1)*A(1,N);
 for K = 2 : M 
 S = S+A(N,K)*A(K,N);
 end;
 A(N,N) = (A(N,N)-S)/XL(N);
 single(A)
 pause 
 if abs(A(N,N)) <= 1.0e-20 
 disp('A = LU but the matrix is singular')
 end
% If A(N,N) = 0 then A = LU but the matrix is singular.
% Process is complete, all entries of A have been determined.
% STEP 7
 L=tril(A,-1) + eye(N); U=triu(A);
 if ISW == 0      
 for K=1:N
 L(K,K)=XL(K);
 end
 else
 for K=1:N
 U(K,K)=XL(K);
 end    
 end

 % MODIFICATIONS
 D = diag(diag(U));
 U = D \ U;

 disp('L='), disp(single(L))
 disp('D='), disp(single(D))
 disp('U='), disp(single(U))
 
 end;
 end;
 
 
 if OK == FALSE 
 fprintf(1,'The matrix does not have an LU factorization.\n');
 end;
 

