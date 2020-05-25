classdef SNPLIB < handle
    %SNPLIB Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Dependent = true)
        nSNPs
    end
    properties
        nThreads = 1
    end
    properties (SetAccess = private)
        nSamples = 0
        GENO
    end
    methods
        function obj = SNPLIB()
            
        end
    end
    methods (Hidden = true)
        function  p = p_mahal_dist(obj,x,y,sigma)
            d = diag([x(:)-1.0,y(:)-1.0]/sigma*[x(:)-1.0,y(:)-1.0]');
            d = d*obj.nSamples/(obj.nSamples-1)^2;
            p = betacdf(d,1,(obj.nSamples-3)/2,'upper');
            p = reshape(p,size(x));
        end
    end
    methods % Gets & Sets
        function out = get.nSNPs(obj)
            out = size(obj.GENO,2);
        end
        function set.nThreads(obj,num_threads)
            if (num_threads<=1)
                obj.nThreads = 1;
            end
            obj.nThreads = num_threads;
        end
    end
    methods
        
    end
    methods % Population Structures
        function af = CalcAlleleFrequency(obj)
            af = CalcAlleleFrequencies_(obj.GENO,obj.nSamples);
        end
        function af = CalcAdjustedAF(obj, scores)
            scores_all = [ones(obj.nSamples,1),scores];
            af = CalcAdjustedAF_(obj.GENO,scores_all,obj.nThreads);
        end
        function maf = CalcAdjustedMAF(obj, scores)
            scores_all = [ones(obj.nSamples,1),scores];
            maf = CalcAdjustedMAF_(obj.GENO,scores_all,obj.nThreads);
        end
        function [matrix,gcta_diag, p_values] = CalcGRMMatrix(obj,check)
            if nargin<2
                check = false;
            end
            af = obj.CalcAlleleFrequency();
            matrix = CalcGRMMatrix_(obj.GENO,af,obj.nSamples,obj.nThreads);
            gcta_diag = CalcGCTADiagonal_(obj.GENO,af,obj.nSamples,obj.nThreads);
            p_values = [];
            if check
                ibs = obj.CalcIBSMatrix();
                x = [diag(matrix)./diag(ibs),gcta_diag./diag(ibs)];
                sigma = robustcov(x);
                % Empirical Mahalanobis Distance follows a beta
                % distribution
                plot(x(:,1),x(:,2),'.');
                hold on
                fun = @(x,y) p_mahal_dist(obj, x,y,sigma);
                x_step = range(x(:,1))/100;
                y_step = range(x(:,2))/100;
                [X,Y] = meshgrid(min(x(:,1)):x_step:max(x(:,1)),min(x(:,2)):y_step:max(x(:,2)));
                Z = fun(X,Y);
                contour(X,Y,Z,[0.001,0.01,0.1],'showtext','on');
                hold off
                p_values = p_mahal_dist(obj,x(:,1),x(:,2),sigma);
            end
        end
        function matrix = CalcGCTAMatrix(obj)
            af = obj.CalcAlleleFrequency();
            matrix = CalcGRMMatrix_(obj.GENO,af,obj.nSamples,obj.nThreads);
            diagonal = CalcGCTADiagonal_(obj.GENO,af,obj.nSamples,obj.nThreads);
            matrix(eye(obj.nSamples)==1) = diagonal;
        end
        function matrix = CalcIBSMatrix(obj)
            matrix = CalcIBSMatrix_(obj.GENO,obj.nSamples,obj.nThreads);
        end
        function matrix = CalcKING(obj)
            matrix = CalcKING_(obj.GENO,obj.nSamples,obj.nThreads);
        end
        function ind = FindUnrelated(obj, threshold)
            if nargin<2
                threshold = 0.044;
            end
            matrix = obj.CalcKING();
            ind = FindUnrelatedGroup_(matrix, threshold);
        end
        function [matrix, gcta_diag, p_values] = CalcAdjustedGRMMatrix(obj,scores,check)
            if nargin < 3
                check = false;
            end
            scores_all = [ones(obj.nSamples,1),scores];
            [matrix, gcta_diag] = CalcAdjustedGRM_(obj.GENO,scores_all,obj.nThreads);
            gcta_diag = diag(matrix) + gcta_diag;
            p_values = [];
            if check
                ibs = obj.CalcIBSMatrix();
                x = [diag(matrix)./diag(ibs),gcta_diag./diag(ibs)];
                sigma = robustcov(x);
                % Empirical Mahalanobis Distance follows a beta
                % distribution
                plot(x(:,1),x(:,2),'.');
                hold on
                fun = @(x,y) p_mahal_dist(obj, x,y,sigma);
                x_step = range(x(:,1))/100;
                y_step = range(x(:,2))/100;
                [X,Y] = meshgrid(min(x(:,1)):x_step:max(x(:,1)),min(x(:,2)):y_step:max(x(:,2)));
                Z = fun(X,Y);
                contour(X,Y,Z,[0.001,0.01,0.1],'showtext','on');
                hold off
                p_values = p_mahal_dist(obj,x(:,1),x(:,2),sigma);
            end
        end
        function [matrix, gcta_diag, p_values] = CalcAdmixedGRMMatrix(obj,pop_af,pop,check)
            if nargin < 4
                check = false;
            end
            [matrix, gcta_diag] = CalcAdmixedGRM_(obj.GENO,pop_af,pop,obj.nThreads);
            gcta_diag = diag(matrix) + gcta_diag;
            p_values = [];
            if check
                ibs = obj.CalcIBSMatrix();
                x = [diag(matrix)./diag(ibs),gcta_diag./diag(ibs)];
                sigma = robustcov(x);
                % Empirical Mahalanobis Distance follows a beta
                % distribution
                plot(x(:,1),x(:,2),'.');
                hold on
                fun = @(x,y) p_mahal_dist(obj, x,y,sigma);
                x_step = range(x(:,1))/100;
                y_step = range(x(:,2))/100;
                [X,Y] = meshgrid(min(x(:,1)):x_step:max(x(:,1)),min(x(:,2)):y_step:max(x(:,2)));
                Z = fun(X,Y);
                contour(X,Y,Z,[0.001,0.01,0.1],'showtext','on');
                hold off
                p_values = p_mahal_dist(obj,x(:,1),x(:,2),sigma);
            end
        end
        function scores = CalcPCAScores(obj,nComponents)
            grm = obj.CalcGRMMatrix();
            [V,D] = eig(grm,'vector');
            [~,ind] = sort(D,'descend');
            scores = V(:,ind(1:nComponents));
        end
        function scores = CalcSUGIBSScores(obj, nComponents)
            ibs = obj.CalcIBSMatrix();
            d = sum(ibs);
            d = diag(d.^(-0.5));
            ugrm = CalcUGRMMatrix_(obj.GENO,obj.nSamples,obj.nThreads);
            I = d*ugrm*d;
            [V,D] = eig(I,'vector');
            [~,ind] = sort(D,'descend');
            scores = d*V(:,ind(2:nComponents+1));
        end
        function loadings = CalcPCALoadingsExact(obj, nComponents)
            af = obj.CalcAlleleFrequency();
            A = UnpackGRMGeno_(obj.GENO,af,obj.nSamples)';
            [U,S,~] = svds(A,nComponents);
            loadings = S\U';
        end
        function loadings = CalcUPCALoadingsExact(obj, nComponents)
            A = UnpackUGeno_(obj.GENO,obj.nSamples)';
            [U,S,~] = svds(A,nComponents);
            loadings = S\U';
        end
        function loadings = CalcSUGIBSLoadingsExact(obj, nComponents)
            ibs = obj.CalcIBSMatrix();
            d = sum(ibs);
            d = diag(d.^(-0.5));
            A = UnpackUGeno_(obj.GENO,obj.nSamples);
            A = A'*d;
            [U,S,~] = svds(A,nComponents+1);
            S = S(2:nComponents+1,2:nComponents+1);
            U = U(:,2:nComponents+1);
            loadings = S\U';
        end
    end
    methods % GWAS
        function [betas,rho2,pvalues] = CalcCCAGWAS(obj, trait)
            Y = trait-repmat(mean(trait,1),[obj.nSamples,1]);
            [betas,rho2] = CalcCCAGWAS_(obj.GENO,Y,obj.nThreads);
            nDims = size(trait,2);
            lambda = 1-rho2;
            t = (nDims^2-4)/(nDims^2+1-5);
            t = sqrt(t);
            w = obj.nSamples-(nDims+4)/2;
            df1 = nDims;
            df2 = w*t-nDims/2+1;
            lambda = lambda.^(1/t);
            F = (1-lambda)./lambda*df2/df1;
            pvalues = fcdf(F,df1,df2,'upper');
        end
        function [rho,pvalues] = CalcCCAReplication(obj,scores,betas)
            rho = CalcCCAReplication_(obj.GENO,scores,betas,obj.nThreads);
            pvalues = tcdf(abs(rho),obj.nSamples-2,'upper')*2;
        end
    end
    methods
        function [SNPs, Samples] = importPLINKDATA(obj,bfile)
            filename = [bfile,'.bim'];
            formatSpec = '%s%s%f%f%s%s';
            fileID = fopen(filename,'r');
            dataArray = textscan(fileID, formatSpec, 'EmptyValue' ,NaN, 'ReturnOnError', false);
            fclose(fileID);
            chr = dataArray{1};
            CHR = str2double(chr);
            CHR(strcmp(chr,'X')) = 23;
            CHR(strcmp(chr,'Y')) = 24;
            CHR(strcmp(chr,'XY')) = 25;
            CHR(strcmp(chr,'MT')) = 26;
            RSID = dataArray{2};
            POS = dataArray{4};
            ALT = upper(dataArray{5});
            REF = upper(dataArray{6});
            SNPs = table(CHR,RSID,POS,ALT,REF);
            filename = [bfile,'.fam'];
            formatSpec = '%s%s%s%s%f%s';
            fileID = fopen(filename,'r');
            dataArray = textscan(fileID, formatSpec, 'MultipleDelimsAsOne', true, 'EmptyValue' ,NaN, 'ReturnOnError', false);
            fclose(fileID);
            FID = dataArray{1};
            IID = dataArray{2};
            PID = dataArray{3};
            MID = dataArray{4};
            Sex = dataArray{5};
            obj.nSamples = length(FID);
            Samples = table(FID,IID,PID,MID,Sex);
            fid = fopen([bfile '.bed'],'rb');
            bin = fread(fid,inf,'uint8=>uint8');
            fclose(fid);
            L = ceil(obj.nSamples/4);
            obj.GENO = reshape(bin(4:end),[L length(REF)]);
        end
        function GenerateIndividuals(obj, af)
            obj.nSamples = size(af,1);
            obj.GENO = GenerateIndividuals_(af);
        end
        function GenerateAdmixedIndividuals(obj, af, num_samples)
            obj.nSamples = num_samples;
            obj.GENO = GenerateAdmixedIndividuals_(af, num_samples);
        end
        function GeneratePairwiseSiblings(obj, parent_obj)
            obj.nSamples = parent_obj.nSamples;
            obj.GENO = GeneratePairwiseSiblings_(parent_obj.GENO,obj.nSamples);
        end
        function extract(obj,ind)
            obj.GENO = obj.GENO(:,ind);
        end
        function FlipGeno(obj,ind)
            FlipGeno_(obj.GENO,obj.nSamples,int32(ind-1));
        end
        function keep(obj,ind)
            Keep_(obj.GENO,obj.nSamples,int32(ind-1));
        end
        function geno_d = UnpackGeno(obj)
            geno_d = UnpackGeno_(obj.GENO,obj.nSamples);
        end
    end
end