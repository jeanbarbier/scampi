classdef FWHTLinTrans_hadamard < LinTrans
    % FWHTLinTrans:  Linear transform class for the fast Walsh-Hadamard
    % Transform
    
    properties
        %Signal model is y = Ax. y and x are vectorizations of a signal
        %which is 1-dimensional
        
        %vector giving the size of the signal x. The signal x is
        %assumed to be the complete output of a FWHT with no samples
        %omitted. Each entry of xSize should be an integer multiple of the
        %corresponding entry of ySize
        xSize;
        
        %vector specifying the m samples of the ySize signal
        %y. These are indices into the 1-dimensional FWHT output with
        %size ySize. These indices can be converted to 1-length subscripts
        %for addressing the array using the ind2sub function.If set to an
        %empty matrix, y is assumed to be the complete ySize FWHT
        %result vectorized.
        ySamples;
        
        %The logarithm of xSize needed for scaling.
        log_xSize;
        
        %actual vector for sign flip
        signVec;

        randP;        
        rp;
        N;            
        M;

    end
    
    methods
        
        % Constructor
        function obj = FWHTLinTrans_hadamard(xSize,ySamples,randP,rp,M)        
            
            %Check for integer multiples
            if mod(log2(xSize),1) ~= 0
                warning('The length of the signal is not a power of 2. The signal will be zero padded to the next highest power of 2.')
            end
            
            %Force ySamples to be a column vector
            ySamples = ySamples(:);
            
            %Set the properties
            obj = obj@LinTrans;
            obj.ySamples = ySamples;
            obj.log_xSize = ceil(log2(xSize));
            obj.xSize = 2^obj.log_xSize;
            obj.randP = randP;
            obj.N = xSize;
            obj.M = M;
            obj.rp = rp;

        end
        
        % Size
        function [m,n] = size(obj)
            
            m = length(obj.ySamples);
            
            %Length of signal x
            n = obj.xSize;
            
            if m>n
                error('The number of measurements cannot be greater than the number of signal elements.')
            end
            
        end
        
        % Replace ySamples with m sample locations chosen uniformly at
        % random
        function ind = ySamplesRandom(obj,m)
            
            %Determine largest possible sample index
            indMax = obj.xSize;
            
            %Generate a random ordering of all possible ySamples indices
            indSamples = randperm(indMax);
            
            %keep only the first m of these random sample indices and
            %ensure that it is a column vector
            ind = reshape(indSamples(1:m),[],1);
        end
        
        %set ySamples using an mx1 matrix of subscripts, rather than
        %providing the indices. This may be more natural for some users.
        function ySamplesSetFromSubScripts(obj,subVals)
            
            %Check size
            if size(subVals,2) ~= 1
                error('Input should be Mx1.')
            end
            
            obj.ySamples = subVals(:);
            
        end
        
        
        % Matrix multiply
        function Z = mult(obj,X)

            numBlockL = 1;
            numBlockC = 1;
            J = 1 / obj.N;
            Nblock = obj.N;
            Mblock = obj.M;
            
            if (isrow(X) ); X = X'; end
            Z = zeros(sum(Mblock), 1);

            % randomization to break the structure
            X(obj.randP) = X; % good

            lastZ = 0;
            for l = 1 : numBlockL
                
                Y = 0;
                for c = 1 : numBlockC
                    if (J(l, c) ~= 0)
                        YY = sign(J(l, c) ) * sqrt(abs(J(l, c) ) ) * hadamards(X((c - 1) * Nblock + 1 : c * Nblock) );
                        Y = Y + YY(obj.rp{l, c}(1 : Mblock(l) ) );
                    end
                end                            
                
                Z(lastZ + 1 : lastZ + Mblock(l) ) = Y;
                lastZ = lastZ + Mblock(l);
            end

            if (isrow(X) ); Z = Z'; end                            
        end
        
        
        % Matrix multiply transpose
        function Z = multTr(obj,X)

            numBlockL = 1;
            numBlockC = 1;
            J = 1 / obj.N;
            Nblock = obj.N;
            Mblock = obj.M;
            
            if (isrow(X) ); X = X'; end

            Z = zeros(numBlockC * Nblock, 1);
            ZZ = Z;
            zero = zeros(Nblock, 1);

            lastZ = 1;
            for c = 1 : numBlockC
                
                u = 1; l_ = 1;
                while (J(l_, c) == 0); u = u + Mblock(l_); l_ = l_ + 1; end
                
                Y = 0; YY = 0;
                for (l = 1 : numBlockL)
                    
                    if (J(l, c) ~= 0)
                        XX = zero;
                        S = X(u : min(end, u + Nblock - 1) );
                        msS = max(size(S) );
                        XX(obj.rp{l, c}(1 : min(msS, Mblock(l) ) ) ) = S(1 : min(msS, Mblock(l) ) );
                        
                        Y = Y + sign(J(l, c) ) * sqrt(abs(J(l, c) ) ) * hadamards(XX);
                        u = u + Mblock(l);                                                
                    end
                end
                
                Z(lastZ : lastZ + Nblock - 1) = Y;                                            
                lastZ = lastZ + Nblock;
            end
           
            if (isrow(X) ); Z = Z'; end

            % derandomization 
            Z = Z(obj.randP); % good            
        end
        
        
        % Matrix multiply with square
        function Z = multSq(obj,X)

            numBlockL = 1;
            numBlockC = 1;
            J = 1 / obj.N;
            Nblock = obj.N;
            Mblock = obj.M;
            
            if (isrow(X) ); X = X'; end
            Z = zeros(sum(Mblock), 1);

            % randomization to break the structure
            X(obj.randP) = X; % good

            lastZ = 0;
            for l = 1 : numBlockL
                
                Y = 0;
                for c = 1 : numBlockC
                    if (J(l, c) ~= 0)
                        YY = abs(J(l, c) ) * sum(X((c - 1) * Nblock + 1 : c * Nblock) ) * ones(Mblock(1), 1);
                        Y = Y + YY(1 : Mblock(l) );
                    end
                end
                
                Z(lastZ + 1 : lastZ + Mblock(l) ) = Y;
                lastZ = lastZ + Mblock(l);
            end

            if (isrow(X) ); Z = Z'; end                    
        end
        
        
        % Matrix multiply transpose
        function Z = multSqTr(obj,X)

            numBlockL = 1;
            numBlockC = 1;
            J = 1 / obj.N;
            Nblock = obj.N;
            Mblock = obj.M;
            
            if (isrow(X) ); X = X'; end
            Z = zeros(numBlockC * Nblock, 1);

            lastZ = 1;
            for c = 1 : numBlockC
                
                u = 1; l_ = 1;
                while (J(l_, c) == 0); u = u + Mblock(l_); l_ = l_ + 1; end
                
                Y = 0;
                for (l = 1 : numBlockL)
                    
                    if (J(l, c) ~= 0)
                        Y = Y + abs(J(l, c) ) * sum(X(u : min(end, u + Mblock(l) - 1) ) ) * ones(Nblock, 1);
                        u = u + Mblock(l);
                    end
                end
                
                Z(lastZ : lastZ + Nblock - 1) = Y;
                lastZ = lastZ + Nblock;
            end

            if (isrow(X) ); Z = Z'; end

            % derandomization 
            Z = Z(obj.randP); % good                
        end                
    end
end
