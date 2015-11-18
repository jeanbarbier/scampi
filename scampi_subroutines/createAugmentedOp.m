function [op, opSq, opT, opSqT, l, L, C, varargout] = createAugmentedOp(hadamard, N, subrate, firstMode)

% inputs :    hadamard = 0/1 : with hadamard ou random measurement matrix? if random, the matrix can be visualized, else the operators are function handle
%             N = n^2 : size of the picture 
%             subrate 
% outputs :   op : direct operator 
%             opSq : squarred op
%             opT : transposed op
%             opT : transposed squarred op
%             L : number of lines of op
%             l : number of lines of the matDiffs
%             C : number of columns of op
%             varargout{1} : random perm (for hadamard op)
%             varargout{2} : flip signs (for hadamard op)
%             varargout{3} : random perm of the signal components (for hadamard op)
%             varargout{4} : the random measurement matrix 
% remark :    all the operators apply on column vector

if N > 64^2; hadamard = 1; end;
M = round(subrate * N);

%% hadamard creation
[rp, fs, rp2] = createRandomLinesAndSignsPermutationForOperators(1, 1, 1 / N, round(subrate * N), N);        
if firstMode; 
    randInd = floor(rand * round(subrate * N) );
    while randInd == 0; randInd = floor(rand * round(subrate * N) ); end
    rp{1, 1}(randInd) = 1; 
end
fs = []; % No sign flip: does not change the performances but faster.        
varargout{1} = rp;
varargout{2} = fs;
varargout{3} = rp2;

%% mappings for messages, nodes are numeroted as im=[1,4,7;2,5,8;3,6,9] and the measured signal is reshape(im,1,[])
for i = 1 : N;
    if i <= N - sqrt(N); map_R(i) = i + sqrt(N);
    else map_R(i) = 0; end
    
    if i > sqrt(N); map_L(i) = i - sqrt(N);
    else map_L(i) = 0; end
    
    if mod(i, sqrt(N) ) ~= 0; map_D(i) = i + 1;
    else map_D(i) = 0; end
    
    if mod(i, sqrt(N) ) ~= 1; map_U(i) = i - 1;
    else map_U(i) = 0; end    
end

%% create the differences matrix
map_1 = map_R;
map_2 = map_D;
RnonZ = map_1 ~= 0;
DnonZ = map_2 ~= 0;

a = [1 : N];
b = [1 : N] + N;
matDiffs = sparse([a(RnonZ), a(RnonZ), b(DnonZ), b(DnonZ)], [map_1(RnonZ), a(RnonZ), map_2(DnonZ), b(DnonZ) - N], [-ones(size(a(RnonZ) ) ), ones(size(a(RnonZ) ) ), -ones(size(b(DnonZ) ) ), ones(size(b(DnonZ) ) ) ], 2 * N, N);

matDiffs = matDiffs(sum(abs(matDiffs.') ) ~= 0, :);
[l, c] = size(matDiffs);
matDiffs = [matDiffs, -speye(l) ];
[l, c] = size(matDiffs); 

%% create the matrix and measure for AMP
if hadamard == 0;
    A = randn(M, N) ./ sqrt(N);
    op = [[A, sparse(M, c - N)]; matDiffs]; 
    opSq = op.^2;
    opT = op.';
    opSqT = opT.^2;   
else  
    A = [];
    matDiffsT = matDiffs.';
    matDiffsSq = matDiffs.^2;
    matDiffsSqT = matDiffsT.^2;  

    op = @(z) [MultSeededHadamard(z(1 : N), 1 / N, 1, 1, M, N, rp, fs, rp2); matDiffs * z];
    opSq = @(z) [mean(z(1 : N) ) * ones(M, 1); matDiffsSq * z];
    opT = @(z) [MultSeededHadamardTranspose(z(1 : M), 1 / N, 1, 1, M, N, rp, fs, rp2) + matDiffsT(1 : N, :) * z(M + 1 : end); matDiffsT(N + 1 : end, :) * z(M + 1 : end) ];    
    opSqT = @(z) [sum(z(1 : M) ) / N + matDiffsSqT(1 : N, :) * z(M + 1 : end); matDiffsSqT(N + 1 : end, :) * z(M + 1 : end) ];
end

varargout{4} = A;
L = M + l;
C = c;

end