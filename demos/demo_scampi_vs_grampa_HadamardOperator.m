% In order to use this demo, the GAMP package (http://sourceforge.net/projects/gampmatlab/) must be installed and in the matlab path 

% set the path to the folder containing the images
pathImages = '/images/';

% choice of the image
N = 512^2; % size of image (power of two with the present Hadamard operator)
image_ = 'lena';

% which algorithms do you want to run?
scampi_ = 1;
grampa_ = 1;

% scampi options
varNoiseDiffs = 1e-1; % initial noise variance associated to the dual variables
opt.learnNoise = 1; % use the Bethe free energy to learn the noise variances
opt.dump_learn = 0.9; % damping for the learning of the noise variance (the larger, the higher the damping)
opt.dump_mes = 0.1; % damping of the AMP algorithm
opt.nb_iter = 300; % max number of iterations
opt.conv = 1e-5; % convergence criterion
opt.print = 10; % print infos every *opt.print* iterations
opt.showDynamics = 1; % show the reconstructed image in real time, every *opt.print* iterations

% compressive imaging recontruction problem definition
subrateTab = [0.1, 0.5]; % subsampling rate(s)
isnrTab = [10, 20, 40, inf]; % input snr(s) (the higher, the less noisy the problem is)
omegaTab = [0 : 5]; % snipe prior parameter(s) (small values are generally better), see http://arxiv.org/abs/1312.3968 for details 
numSamples = 1; % number of samples (i.e. different instances of measurement matrix and noise) reconstructed for each settings
mismatchTab = 10.^[-4 : 4]; % values of noise mismatch, i.e. how wrong is the initial noise value given to the algorithm(s)

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

% main loop
for i = 1 : numel(isnrTab)
    isnr = isnrTab(i);

    for s = 1 : numel(subrateTab)
        subrate = subrateTab(s);
        
        for o = 1 : numel(omegaTab)
            opt.omega = omegaTab(o);            

            for m = 1 : numel(mismatchTab)
                mismatch = mismatchTab(m);

                for t = 1 : numSamples

                    % create the augmented operators for the cosparse TV analysis model with dual variables
                    [op, opSq, opTr, opSqTr, l, ~, Ntot, rp, ~, rp2] = createAugmentedOp(N, subrate, 1);

                    % compressive linear measurement
                    M = round(subrate * N);
                    y = MultSeededHadamard(opt.signal, 1 / N, 1, 1, M, N, rp, [], rp2);                

                    % noise
                    nvar = mean(y(y < max(abs(y) ) ).^2) * exp(-log(10) * isnr / 10); % noise variance (avoiding to take into account the measurement associated to the first Hadamard mode)
                    noise = sqrt(nvar) .* randn(round(subrate * N), 1); % zero mean additive white Gaussian noise    
                    y = y + noise; 
                    yz = [y; sparse(l, 1) ]; % augmented measurement               

                    % scampi
                    if scampi_                             
                        if mismatch ~= 0; nvar_sc = mismatch * nvar; else nvar_sc = nvar; end
                        if opt.learnNoise; % noise variance (associated to pixels and dual variables)    
                            opt.var_noise = [1e-1 * ones(M, 1); varNoiseDiffs * ones(l, 1) ];
                        else 
                            opt.var_noise = [nvar_sc * ones(M, 1); varNoiseDiffs * ones(l, 1) ];
                        end                         
                        opt.M = numel(yz); % size of augmented measurement
                        opt.N = Ntot; % size of augmented signal
                        opt.Ntrue = N; % size of image
                        opt.Mtrue = M; % size of original measurement
                                                                                                                    
                        tic; X_scampi = scampi_solver(yz, op, opTr, opSq, opSqTr, opt); time_scampi = toc;
                        nsnr_scampi = 10 * log10(norm(opt.signal(:) - X_scampi(1 : N) ).^2 ./ norm(opt.signal).^2);                
                        mse_scampi = mean((X_scampi(1 : N) - opt.signal(:) ).^2);                

                        results_scampi_mse{i, s, o, m, t} = mse_scampi;
                        results_scampi_nsnr{i, s, o, m, t} = nsnr_scampi;
                        results_scampi_time{i, s, o, m, t} = time_scampi;
                    end

                    % grampa
                    if grampa_    
                        RandProj = FWHTLinTrans_hadamard(N, randperm(N, M), rp2, rp, M);
                        TVProj = LinTransTV([sqrt(N), sqrt(N) ], 1);
                        tmp{1} = RandProj;
                        tmp{2} = TVProj;
                        AnalysisProj = LinTransConcat(tmp');
                        GrampaOptions = GampOpt('legacyOut',false,'uniformVariance',true,'adaptStep',false);
                        GrampaOptions.xvar0 = 1;
                        GrampaOptions.step = .1; % initial stepsize
                        GrampaOptions.stepMax = .5; % maximum stepsize
                        GrampaOptions.stepIncr = 1.01; % stepsize increase rate
                        GrampaOptions.tol = opt.conv; % stopping tolerance
                        GrampaOptions.nit = opt.nb_iter; % maximum number of iterations
                        GrampaOptions.varNorm = true; % turn on internal normalization
                        GrampaOptions.zvarToPvarMax = inf; % analysis variance clipping ratio
                        GrampaOptions.adaptStepBethe = 1; % use adaptative learning using the Bethe free energy
                        GrampaOptions.verbose = 1;                        
                        GrampaEstimIn = NullEstimIn(0,1);
                        GrampaLinTrans = AnalysisProj;
                        if mismatch ~= 0; nvar_gr = mismatch * nvar; else nvar_gr = nvar; end
                        MeasEstimOut = AwgnEstimOut(y,nvar_gr);    
                        AnaEstimOut = SNIPEstim(opt.omega);
                        GrampaEstimOut = EstimOutConcat({MeasEstimOut;AnaEstimOut},[M,size(TVProj,1)]);
                        tic; estFin1 = gampEst(GrampaEstimIn,GrampaEstimOut,GrampaLinTrans,GrampaOptions); time_grampa = toc;
                        X_grampa = estFin1.xhat;                        
                        nsnr_grampa = 10 * log10(norm(opt.signal(:) - X_grampa(:) ).^2 ./ norm(opt.signal).^2);
                        mse_grampa = mean((X_grampa(:) - opt.signal(:) ).^2);  

                        results_grampa_mse{i, s, o, m, t} = mse_grampa;
                        results_grampa_nsnr{i, s, o, m, t} = nsnr_grampa;
                        results_grampa_time{i, s, o, m, t} = time_grampa;
                    end  

                    results_nvar{i, s, o, m, t} = nvar;                
                    disp('isnr  subrate   mismatch     nsnr_scampi     nsnr_grampa'); pr = sprintf('%d   %0.1f       %0.1e      %0.2e         %0.2e', [isnr, subrate, mismatch, nsnr_scampi, nsnr_grampa] ); disp(pr);            
                end

            end

        end

    end

end

% plot the charts comparing scampi and grampa
plotFigs;