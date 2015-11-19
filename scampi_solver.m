function [X_hat, Diffs_hat] = scampi_solver(Y, op, opT, opSq, opSqT, opt)

% initialisation 
Y = Y(:); 
opt.signal = opt.signal(:);

W = zeros(opt.M, 1); 
W_new = W;
R_init = zeros(opt.N, 1); 
S2_init = R_init;
S2_init(1 : opt.Ntrue) = opt.var_noise(1); 
S2_init(1 + opt.Ntrue : end) = opt.var_noise(end); 
av_mess_init = R_init;
var_mess_init = S2_init; 
V = mean(var_mess_init) / (opt.Mtrue / opt.Ntrue) * ones(opt.M, 1); 
V_new = V; 

prior = PriorScampi(opt.Ntrue, R_init, S2_init, av_mess_init, var_mess_init, opt.omega);

% Starting main code
t = 1;
while t <= opt.nb_iter

    av_mess_temp = prior.av_mess;
    var_mess_temp = prior.var_mess;    
    
    % AMP for the augmented system
    V_new = opSq(prior.var_mess);    
    W_new = dumping(W, op(prior.av_mess) - ((Y - W) ./ (opt.var_noise + V) ) .* V_new, opt.dump_mes);
    V_new = dumping(V, V_new, opt.dump_mes);   
    var_1 = opSqT(1 ./ (opt.var_noise + V_new) );
    var_2 = opT((Y - W_new) ./ (opt.var_noise + V_new) ); 
    prior.R = var_2 ./ var_1 + prior.av_mess; 
    prior.S2 = 1 ./ var_1;     
    prior = Prior(prior); 
    V = V_new; 
    W = W_new;    
    
    % use the Bethe free energy to learn the noise variance terms    
    if opt.learnNoise;       
        par = (Y - op(prior.av_mess) ).^2;        
        newNoise = .5 .* (par + sqrt(par.^2 + 4 .* opSq(prior.var_mess) .* par) );
        newNoise1 = opt.dump_learn * opt.var_noise(1) + (1 - opt.dump_learn) * mean(newNoise(1 : opt.Mtrue) );    
        newNoise2 = opt.dump_learn * opt.var_noise(end) + (1 - opt.dump_learn) * mean(newNoise(1 + opt.Mtrue : end) );                
        opt.var_noise(1 : opt.Mtrue) = newNoise1;
        opt.var_noise(1 + opt.Mtrue : end) = newNoise2;    
    end     

    % mse and convergence        
    if mod(t, opt.print) == 0
        conv_ = mean(abs(av_mess_temp - prior.av_mess) );    
        if (t > 0) & (conv_ < opt.conv)
            disp('converged'); 
            break; 
        end          
        nsnr_ = 10 * log10(norm(opt.signal(:) - prior.av_mess(1 : opt.Ntrue) ).^2 ./ norm(opt.signal).^2);       
        disp(sprintf('t=%d, convergence=%0.2e, pixel noise var=%0.2e, dual noise var=%0.2e, nsnr=%0.2e', [t, conv_, opt.var_noise(1), opt.var_noise(end), nsnr_] ) );         
    end   

    % show real time reconstruction
    if opt.showDynamics && (mod(t, opt.print) == 0);
        figure(1); 
        imshow(reshape(prior.av_mess(1 : opt.Ntrue), sqrt(opt.Ntrue), sqrt(opt.Ntrue) ) ); 
        drawnow;
    end

    t = t + 1;    
end

X_hat = prior.av_mess(1 : opt.Ntrue);
Diffs_hat = prior.av_mess(opt.Ntrue + 1 : end);

end