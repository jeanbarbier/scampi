function [op, opSq, opT, opSqT, l, L, C, varargout] = createAugmentedOp(N, subrate, firstMode, varargin)

% Two possibilities for the measurement matrix used for the augmented operator : 1) a matrix F is given, 2) if not, a subsampled Hadamard fast operator is used (default).
%
% inputs :    N : total number of pixels of the picture (of size N = n * n)
%             subrate : subsampling rate <= 1
%             firstMode : enforce the precense of the all 1's mode in the Hadamard operator, if used
%             varargin{1} = F : the original measurement matrix             
% outputs :   op : direct operator 
%             opSq : squarred op
%             opT : transposed op
%             opSqT : transposed squarred op
%             l : number of lines of the matDiffs            
%             L : number of lines of op             
%             C : number of columns of op
%             varargout{1} : random perm (for Hadamard op)
%             varargout{2} : fliped signs (for Hadamard op)
%             varargout{3} : random perm of the signal components (for Hadamard op)
% remark :    all the operators apply on column vector

M = round(subrate * N);

if nargin == 3
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
end

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
matDiffsT = matDiffs.';
matDiffsSq = matDiffs.^2;
matDiffsSqT = matDiffsT.^2;  

if nargin == 3
    op = @(z) [MultSeededHadamard(z(1 : N), 1 / N, 1, 1, M, N, rp, fs, rp2); matDiffs * z];
    opSq = @(z) [mean(z(1 : N) ) * ones(M, 1); matDiffsSq * z];
    opT = @(z) [MultSeededHadamardTranspose(z(1 : M), 1 / N, 1, 1, M, N, rp, fs, rp2) + matDiffsT(1 : N, :) * z(M + 1 : end); matDiffsT(N + 1 : end, :) * z(M + 1 : end) ];    
    opSqT = @(z) [sum(z(1 : M) ) / N + matDiffsSqT(1 : N, :) * z(M + 1 : end); matDiffsSqT(N + 1 : end, :) * z(M + 1 : end) ];
else
    op_ = [[varargin{1}, sparse(M, c - N)]; matDiffs]; 
    op = @(z) op_ * z(:);

    opSq_ = op_.^2;
    opSq = @(z) opSq_ * z(:);

    opT_ = op_.';
    opT = @(z) opT_ * z(:);

    opSqT_ = opT_.^2; 
    opSqT = @(z) opSqT_ * z(:);
end 

L = M + l;
C = c;

end