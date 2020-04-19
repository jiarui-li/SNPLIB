classdef SNPLIB < handle
    %SNPLIB Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Dependent = true)
        nSNPs
        nSamples
    end
    properties
        nThreads = 1
    end
    properties (SetAccess = private)
        CHR
        RSID
        POS
        REF
        ALT
        IID
        FID
        PID
        MID
        Sex
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
            out = length(obj.RSID);
        end
        function out = get.nSamples(obj)
            out = length(obj.FID);
        end
        function set.nThreads(obj,num_threads)
            if (num_threads<=1)
                obj.nThreads = 1;
            end
            obj.nThreads = num_threads;
        end
    end
    methods
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
            
            matrix(eye(obj.nSamples)==1) = diagonal;
        end
        function matrix = CalcIBSMatrix(obj)
            matrix = CalcIBSMatrix_(obj.GENO,obj.nSamples,obj.nThreads);
        end
        function matrix = CalcKING(obj)
            matrix = CalcKING_(obj.GENO,obj.nSamples,obj.nThreads);
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
            A = UnpackPCA_(obj.GENO,af,obj.nSamples)';
            [U,S,~] = svds(A,nComponents);
            loadings = S\U';
        end
        function loadings = CalcUPCALoadingsExact(obj, nComponents)
            A = UnpackUG_(obj.GENO,obj.nSamples)';
            [U,S,~] = svds(A,nComponents);
            loadings = S\U';
        end
        function loadings = CalcSUGIBSLoadingsExact(obj, nComponents)
            ibs = obj.CalcIBSMatrix();
            d = sum(ibs);
            d = diag(d.^(-0.5));
            A = UnpackUG_(obj.GENO,obj.nSamples);
            A = A'*d;
            [U,S,~] = svds(A,nComponents+1);
            S = S(2:nComponents+1,2:nComponents+1);
            U = U(:,2:nComponents+1);
            loadings = S\U';
        end
    end
    methods
        function createSIM(obj,geno,N)
            obj.GENO = geno;
            V = size(geno,2);
            obj.CHR = ones(V,1);
            obj.RSID = (1:V)';
            obj.POS = (1:V)';
            obj.REF = repmat({'A'},[V,1]);
            obj.ALT = repmat({'B'},[V,1]);
            obj.FID = (1:N)';
            obj.IID = (1:N)';
            obj.PID = zeros(N,1);
            obj.MID = zeros(N,1);
            obj.Sex = ones(N,1);
        end
        function geno = UnpackGENO(obj)
            geno = UnpackGENO_(obj.GENO,obj.nSamples);
        end
        function extract(obj,ind)
            obj.GENO = obj.GENO(:,ind);
            obj.RSID = obj.RSID(ind);
            obj.CHR = obj.CHR(ind);
            obj.POS = obj.POS(ind);
            obj.REF = obj.REF(ind);
            obj.ALT = obj.ALT(ind);
        end
        function exportPLINKDATA(obj, bfile)
            % EXPORTPLINKDATA Export to Plink binary data
            %   EXPORTPLINKDATA(obj, bfile) export to plink binary data
            %   (bfile.bed, bfile.bim, bfile.fam) from obj
            T = table(obj.CHR,obj.RSID,zeros(obj.nSNPs,1),obj.POS,obj.ALT, obj.REF);
            writetable(T, [bfile,'.txt'], 'WriteVariableNames', false, 'Delimiter', '\t');
            movefile([bfile,'.txt'],[bfile,'.bim']);
            T = table(obj.FID,obj.IID,obj.PID,obj.MID, obj.Sex, zeros(obj.nSamples,1));
            writetable(T, [bfile,'.txt'], 'WriteVariableNames', false, 'Delimiter', '\t');
            movefile([bfile,'.txt'],[bfile,'.fam']);
            fid = fopen([bfile '.bed'],'wb');
            fwrite(fid,[108, 27,1],'uint8');
            fwrite(fid,obj.GENO,'uint8');
            fclose(fid);
        end
        %         function AlignSNPs(obj,ref_obj)
        %             if nargin<2 || ~isa(ref_obj,'SNPLIB')
        %                 error('Please provide a SNPLIB class object as reference!');
        %             end
        %             chrs = unique(ref_obj.CHR);
        %             ind = [];
        %             alleles = union(union(cellstr(obj.ALT),cellstr(ref_obj.ALT)),union(cellstr(obj.REF),cellstr(ref_obj.REF)));
        %             AMAP = containers.Map(alleles,1:length(alleles));
        %             for i=1:length(chrs)
        %                 ind1 = find(obj.CHR==chrs(i));
        %                 ind2 = find(ref_obj.CHR==chrs(i));
        %                 POS1 = uint64(obj.POS(ind1));
        %                 POS2 = uint64(ref_obj.POS(ind2));
        %                 set1 = POS1+...
        %                     uint64(bitshift(cell2mat(AMAP.values(obj.ALT(ind1))),32,'uint64'))+...
        %                     uint64(bitshift(cell2mat(AMAP.values(obj.REF(ind1))),48,'uint64'));
        %                 set1_r = POS1+...
        %                     uint64(bitshift(cell2mat(AMAP.values(obj.REF(ind1))),32,'uint64'))+...
        %                     uint64(bitshift(cell2mat(AMAP.values(obj.ALT(ind1))),48,'uint64'));
        %                 set2 = POS2+...
        %                     uint64(bitshift(cell2mat(AMAP.values(ref_obj.ALT(ind2))),32,'uint64'))+...
        %                     uint64(bitshift(cell2mat(AMAP.values(ref_obj.REF(ind2))),48,'uint64'));
        %                 ind_1 = find(ismember(set1,set2));
        %                 ind_2 = find(ismember(set1_r,set2));
        %                 ind = [ind;ind1(union(ind_1,ind_2))];
        %             end
        %             obj.extract(ind);
        %         end
        %         function AlignAlleles(obj,ref_obj)
        %             % ALIGN_ALLELE Align the alleles with reference dataset
        %             %   ALIGN_ALLELE(obj, ref_obj)
        %             if nargin<2 || ~isa(ref_obj,'SNPLIB')
        %                 error('Please provide a SNPLIB class object as reference!');
        %             end
        %             ind = find(~strcmp(obj.REF,ref_obj.REF));
        %             FlipGeno_(obj.GENO,obj.nSamples,int32(ind-1));
        %         end
        function  Align(obj,ref_obj)
            % Assuming two datasets have the same SNP set
            ind = find(~strcmp(obj.ALT,ref_obj.ALT));
            obj.FlipGeno(ind);
            obj.RSID = ref_obj.RSID;
            obj.ALT = ref_obj.ALT;
            obj.REF = ref_obj.REF;
        end
        function FlipGeno(obj,ind)
            FlipGeno_(obj.GENO,obj.nSamples,int32(ind-1));
        end
        function keep(obj,ind)
            obj.FID = obj.FID(ind);
            obj.IID = obj.IID(ind);
            obj.PID = obj.PID(ind);
            obj.MID = obj.MID(ind);
            obj.Sex = obj.Sex(ind);
            Keep_(obj.GENO,obj.nSamples,int32(ind-1));
        end
    end
end