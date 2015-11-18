% set the path to the folder containing the images
pathImages = '~/Desktop/scampi/images/';

% choice of the image
N = 512^2; % size of image (power of two with Hadamard operator)
image_ = 'lena';

% compressive imaging recontruction problem definition
subrate = 0.1; % subsampling rate
isnr = 20; % input snr (the higher, the less noisy the problem is)
opt.omega = 0; % snipe prior parameter (small values are generally better), see http://arxiv.org/abs/1312.3968 for details 

% scampi options
varNoiseDiffs = 1e-1; % initial noise variance associated to the dual variables
opt.learnNoise = 1; % use the Bethe free energy to learn the noise variances
opt.dump_learn = 0.9; % damping for the learning of the noise variance (the larger, the higher the damping)
opt.dump_mes = 0.1; % damping of the AMP algorithm
opt.nb_iter = 300; % max number of iterations
opt.conv = 1e-5; % convergence criterion
opt.print = 10; % print infos every *opt.print* iterations
opt.showDynamics = 0; % show the reconstructed image in real time

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

% image
switch image_
    case 'cameraman'
        original_image = double(imread([pathImages, 'cam512.jpg'] ) );   
    case 'boat'
        original_image = double(imread([pathImages, 'boat512.png'] ) ); 
    case 'lena'
        original_image = double(imread([pathImages, 'lena512.bmp'] ) );
    case 'barbara'
        original_image = double(imread([pathImages, 'barbara512.png'] ) );       
    case 'peppers'
        original_image = double(imread([pathImages, 'peppers512.jpg'] ) ); 
    case 'bridge'
        original_image = double(imread([pathImages, 'bridge512.jpg'] ) );
    case 'mandrill'
        original_image = double(imread([pathImages, 'mandrill512.jpg'] ) );  
    case 'phantom'    
        original_image = phantom(sqrt(N) );    
    otherwise           
        disp('unknown image');
end
downscale = sqrt(N) / max(size(original_image) );
original_image = round(imresize(original_image, round(downscale * size(original_image) ) ) ); % downscale image if needed  
original_image = rescaleImage(original_image); % rescale image to have pixel values in [0,1]
opt.signal = original_image(:);

% create the augmented operators for the cosparse TV analysis model with dual variables
[op, opSq, opTr, opSqTr, l, ~, Ntot, rp, ~, rp2, ~] = createAugmentedOp(1, N, subrate, 1);

% compressive linear measurement
M = round(subrate * N);
y = MultSeededHadamard(opt.signal, 1 / N, 1, 1, M, N, rp, [], rp2);                

% noise
nvar = mean(y(y < max(abs(y) ) ).^2) * exp(-log(10) * isnr / 10); % noise variance (avoiding to take into account the measurement associated to the first Hadamard mode)
noise = sqrt(nvar) .* randn(round(subrate * N), 1); % zero mean additive white Gaussian noise    
y = y + noise; 
yz = [y; sparse(l, 1) ]; % augmented measurement
                             
% some fields required by scampi              
if opt.learnNoise; % noise variance (associated to pixels and dual variables)    
    opt.var_noise = [1e-1 * ones(M, 1); varNoiseDiffs * ones(l, 1) ];
else 
    opt.var_noise = [nvar * ones(M, 1); varNoiseDiffs * ones(l, 1) ];
end                
opt.M = numel(yz); % size of augmented measurement
opt.N = Ntot; % size of augmented signal
opt.Ntrue = N; % size of image
opt.Mtrue = M; % size of original measurement
       
% scampi                                                                                                             
[X_scampi, D_scampi] = scampi_solver(yz, op, opTr, opSq, opSqTr, opt); % X_scampi is the reconstructed image, D_scampi, the estimation of the dual variables
nsnr = 10 * log10(norm(opt.signal(:) - X_scampi(1 : N) ).^2 ./ norm(opt.signal).^2); % final nsnr
disp(sprintf('final nsnr=%0.2e', nsnr) );
subplot(1, 2, 1); imshow(original_image); title('original', 'Fontsize', 15);
subplot(1, 2, 2); imshow(reshape(X_scampi, sqrt(N), sqrt(N) ) ); title('scampi estimate', 'Fontsize', 15);