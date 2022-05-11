function [fMatrix] = FfromEightPnts(xs, xss)

% xs = p1;
% xss = p2;

%Check points are in correct format
[c1, n1]=size(xs);
[c2, n2]=size(xss);

if ((c1 ~=3) || (c2 ~=3))
    error('Points in incorrect format');
end

if ((n1 < 8) || (n2 < 8))
    error('Not Enough Points');
end

% Normalize the points
Norm_Matrix1 = GetNormMat(xs);
Norm_Matrix2 = GetNormMat(xss);
% normalize the point set.
NormlizedPts1 = Norm_Matrix1 * xs;
NormlizedPts2 = Norm_Matrix2 * xss;

% % Without Normalization
% NormlizedPts1 = xs;
% NormlizedPts2 = xss;


% Compute the essential?? matrix E from point correspendences
% (NormlizedPts1.') * F_q * NormlizedPts2 = 0
% Ax=0

% Cosntruct Matrix A by Kronecker Product
for n =1 : size(xs,2) % 27 columns
    A(n,:) = kron(NormlizedPts2(:,n).',NormlizedPts1(:,n).');
%     B(n,:) = kron(xss(n,:),xs(n,:));
%     a = [1 2 3];
%     b = [4 5 6];
%     kron(a,b)
end

% Find the solution of f by SVD. Af =0 
[U_A,S_A,V_A] = svd(A);
f = V_A(:,end);
F =  reshape(f,3,3);
F = F.';

% Enforce F to Rank 2
[U_F,S_F,V_F] = svd(F);
S_F(3,3) = 0;
F = U_F * S_F * V_F.';

% % Denormalization
fMatrix = Norm_Matrix2.' * F * Norm_Matrix1;
% fMatrix =  F ;

end


