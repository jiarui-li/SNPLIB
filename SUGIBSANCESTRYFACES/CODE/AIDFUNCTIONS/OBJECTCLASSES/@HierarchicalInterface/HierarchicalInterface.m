classdef HierarchicalInterface < superClassLight
    properties
        nL = 6;
    end
    properties (Dependent = true)
        nLC;
    end
    properties (Hidden = true)
        Levels;
        Clusters;
    end
    properties(Hidden = true, Dependent = true)
        indLabels;
        lcLabels;
    end
    methods % CONSTRUCTOR
        function obj = HierarchicalInterface(varargin)
            %obj = obj@superClassLight(varargin{:});
            getLevels(obj);
            getClusters(obj);        
        end
    end
    methods % GENERAL GETTING
        function out = get.nLC(obj)
                 out = 2^obj.nL-1; 
        end
        function out = get.indLabels(obj)
                 out = cell(1,obj.nLC);
                 for i=1:1:obj.nLC
                    out{i} = num2str(i);
                 end
        end
        function out = get.lcLabels(obj)
                 out = cell(1,obj.nLC);
                 for i=1:1:obj.nLC
                    [l,c] = Ind2LC(obj,i);
                    out{i} = [num2str(l) '/' num2str(c)];
                 end
        end
    end
    methods % GENERAL SETTING
       function obj = set.nL(obj,in)
           obj.nL = single(in);
           getLevels(obj);
           getClusters(obj);
       end
    end
    methods % BASIC INTERFACING
       function out = getLevels(obj)
            out = zeros(1,obj.nLC);
            counter = 1;
            for i=1:1:obj.nL
               cl = i*ones(1,2^(i-1));
               out(counter:counter+length(cl)-1) = cl;
               counter = counter+length(cl);
            end
            %out = HierarchicalInterface.convertUInt(out);
            obj.Levels = single(out);
       end
       function out =getClusters(obj)
           out = zeros(1,obj.nLC);
           counter = 1;
           for i=1:1:obj.nL
               cl = 1:2^(i-1);
               out(counter:counter+length(cl)-1) = cl;
               counter = counter+length(cl);
           end
           %out = HierarchicalInterface.convertUInt(out);
           obj.Clusters = single(out);
       end
       function i = LC2Ind(obj,l,c)
           if ~isscalar(l)% multiple conversions to make
              n = length(l);
              i = zeros(1,n);
              for k=1:1:n
                 i(k) = LC2Ind(obj,l(k),c(k)); 
              end
              return; 
           end
           if l > obj.nL; i = []; return; end
           if c > 2^(l-1), i = []; return; end
           i = sum(2.^((1:1:l-1)-1));
           i = i + c;
       end
       function [l,c] = Ind2LC(obj,i)
           if ~isscalar(i)% multiple conversions to make
              n = length(i);
              l = zeros(1,n);
              c = zeros(1,n);
              for k=1:1:n
                 [l(k),c(k)] = Ind2LC(obj,i(k)); 
              end
              return; 
           end
           l = obj.Levels(i);
           c = obj.Clusters(i); 
       end
       function [pi,pl,pc] = getParent(obj,l,c)
          if nargin<3,[l,c] = Ind2LC(obj,l);end%Index was given
          if l==1, pl = [];pc = []; pi= []; return; end
          if c > 2^(l-1), pl = [];pc = []; pi= []; return; end
          pl = l-1;
          pc = ceil(c/2);
          pi = LC2Ind(obj,pl,pc);
       end
       function [ci,cl,cc] = getChildren(obj,l,c)
          if nargin<3,[l,c] = Ind2LC(obj,l);end%Index was given
          if l >= obj.nL; cl=[];cc=[];ci=[];return; end
          if c > 2^(l-1), cl=[];cc=[];ci=[];return; end
          cl = (l+1)*ones(1,2);
          cc = [(2*c)-1 2*c];
          ci = zeros(1,2);
          for k=1:1:2
              ci(k) = LC2Ind(obj,cl(k),cc(k));
          end
       end
       function [pi,pl,pc] = getAllParents(obj,l,c)
          if nargin<3,[l,c] = Ind2LC(obj,l);end%Index was given
          if l==1, pl = [];pc = []; pi= []; return; end
          nrParents = l-1;
          pi = zeros(1,nrParents);pl = zeros(1,nrParents);pc = zeros(1,nrParents);
          currentl = l;currentc = c;
          for k=1:1:nrParents
              [pi(k),pl(k),pc(k)] = getParent(obj,currentl,currentc);
              currentl = pl(k);currentc = pc(k);
          end
       end
       function [ci,cl,cc] = getAllChildren(obj,l,c)
          if nargin<3,[l,c] = Ind2LC(obj,l);end%Index was given 
          if l >= obj.nL; cl=[];cc=[];ci=[];return; end
          if c > 2^(l-1), cl=[];cc=[];ci=[];return; end 
          [di,dl,dc] = getChildren(obj,l,c);
          [d1i,d1l,d1c] = getAllChildren(obj,dl(1),dc(1));
          [d2i,d2l,d2c] = getAllChildren(obj,dl(2),dc(2));
          cl = [dl d1l d2l];
          cc = [dc d1c d2c];
          ci = [di d1i d2i];
       end
       function out = List2Matrix(obj,in)
            nCl = 2^(obj.nL-1); 
            out = zeros(obj.nL,nCl);
            for c=1:nCl
                %c=1;
                [~,~,pc] = getAllParents(obj,obj.nL,c);
                parentEvolution = flip([c pc]);
                for l=1:obj.nL
                    %l=1;
                    ind = LC2Ind(obj,l,parentEvolution(l));
                    out(l,c) = in(ind);
                end
            end
       end
    end
    methods % ADVANCED INTERFACING
       function out = imageMatrix(obj,in,plotFrame)
                if nargin<3, plotFrame = true; end
                out = figure;imagesc(List2Matrix(obj,in));
                nClusters = 2^(obj.nL-1);
                if plotFrame == true,
                    hold on
                    plot([0.5 nClusters+0.5],[0.5:obj.nL+0.5; 0.5:obj.nL+0.5]','-k','LineWidth',2)
                    for level = 1:obj.nL,
                        plot([0.5: nClusters/(2^(level-1)) :nClusters+0.5 ; 0.5: nClusters/(2^(level-1)) :nClusters+0.5],[(level-1) + 0.5 , level + 0.5],'-k','LineWidth',3)
                    end
                end
       end
       function out = getCellData(obj,in,modind,fName)
           if ~iscell(in), out = in;return;end
           if nargin<4, fName = []; end
           % RECURSIVE PULLING OUT OF DATA
           if ~isscalar(modind)
              out = [];
              for i=1:1:length(modind)
                 out = [out, getCellData(obj,in,modind(i))];  %#ok<AGROW>
              end
              return;
           end
           % SINGLE PULL OF DATA
           % determining field name
           if isempty(fName),fNames = fieldnames(in{1});fName = fNames{1};end
           if isscalar(fName),fNames = fieldnames(in{1});fName = fNames{fName};end
           [modL,modC] = Ind2LC(obj,modind);
           eval(['out = in{' num2str(modL) '}.' fName '{' num2str(modC) '};']);
       end
       function out = reduceCellData(obj,in,keep,fName)
           if ~iscell(in), out = in(index,:);return;end
           if nargin<4, fName = []; end
           out = in;
           % determining field name
           if isempty(fName),fNames = fieldnames(in{1});fName = fNames{1};end
           if isscalar(fName),fNames = fieldnames(in{1});fName = fNames{fName};end
           for i = 1:1:obj.nLC
               %i=1;
               [l,c] = Ind2LC(obj,i);
               eval(['tmp = in{l}.' fName '{c};']);
               eval(['out{l}.' fName '{c} = tmp(keep,:);'])
           end
       end
       function out = mergeCellData(obj,in1,in2,fName)
           if ~iscell(in1), out = [in1;in2];return;end
           if nargin<4, fName = []; end
           out = in1;
           % determining field name
           if isempty(fName),fNames = fieldnames(in1{1});fName = fNames{1};end
           if isscalar(fName),fNames = fieldnames(in1{1});fName = fNames{fName};end
           for i = 1:1:obj.nLC
               %i=1;
               [l,c] = Ind2LC(obj,i);
               eval(['tmp1 = in1{l}.' fName '{c};']);
               eval(['tmp2 = in2{l}.' fName '{c};']);
               eval(['out{l}.' fName '{c} = [tmp1;tmp2];'])
           end
       end
       function [out1,out2] = getCellDataSize(obj,in,modind)
           if ~iscell(in), [out1,out2] = size(in);return;end
           % RECURSIVE PULLING OUT OF DATA
           if ~isscalar(modind)
              out1 = zeros(1,length(modind));
              out2 = zeros(1,length(modind));
              for i=1:1:length(modind)
                 [out1(i),out2(i)] = getCellDataSize(obj,in,modind(i)); 
              end
              return;
           end
           % SINGLE PULL OF DATA
           [out1,out2] = size(getCellData(obj,in,modind));
       end
       function out = expandWithChildren(obj,index)
           out = [];
           for i=length(index):-1:1
               %i=length(index);
               children = [index(i) getAllChildren(obj,index(i))];
               out = union(out,children);
           end
       end
       function out = expandWithParents(obj,index)
           out = [];
           for i=length(index):-1:1
               %i=length(index);
               parents = [index(i) getAllParents(obj,index(i))];
               out = union(out,parents);
           end
           
       end
    end
    methods (Static = true)
        function out = convertUInt(in)
           if isempty(in), out = []; return; end
           M = max(in(:));
           if M<=intmax('uint8'),out = uint8(in);return;end
           if M<=intmax('uint16'), out = uint16(in);return;end
           if M<=intmax('uint32'), out = uint32(in);return;end
           out = uint64(in);
       end
    end
end


% function [pl,pc,pi] = getParent(obj,l,c)
%           if nargin<3,[l,c] = Ind2LC(obj,l);end%Index was given
%           if l==1, pl = [];pc = []; pi= []; return; end
%           if c > 2^(l-1), pl = [];pc = []; pi= []; return; end
%           pl = l-1;
%           pc = ceil(c/2); 
%           if nargout<3, return; end
%           pi = LC2Ind(obj,pl,pc);
%        end
%        function [cl,cc,ci] = getChildren(obj,l,c)
%           if nargin<3,[l,c] = Ind2LC(obj,l);end%Index was given
%           if l >= obj.nL; cl=[];cc=[];ci=[];return; end
%           if c > 2^(l-1), cl=[];cc=[];ci=[];return; end
%           cl = (l+1)*ones(1,2);
%           cc = [(2*c)-1 2*c];
%           if nargout<3, return; end
%           ci = zeros(1,2);
%           for k=1:1:2
%               ci(k) = LC2Ind(obj,cl(k),cc(k));
%           end
%        end
%        function [pl,pc,pi] = getAllParents(obj,l,c)
%           if nargin<3,[l,c] = Ind2LC(obj,l);end%Index was given
%           if l==1, pl = [];pc = []; pi= []; return; end
%           nrParents = l-1;
%           pl = zeros(1,nrParents);pc = zeros(1,nrParents);
%           currentl = l;currentc = c;
%           for k=1:1:nrParents
%               [pl(k),pc(k)] = getParent(obj,currentl,currentc);
%               currentl = pl(k);currentc = pc(k);
%           end
%           if nargout<3, return; end
%           pi = zeros(1,nrParents);
%           for k=1:1:nrParents
%              pi(k) =  LC2Ind(obj,pl(k),pc(k));
%           end
%        end
%        function [cl,cc,ci] = getAllChildren(obj,l,c)
%           if nargin<3,[l,c] = Ind2LC(obj,l);end%Index was given 
%           if l >= obj.nL; cl=[];cc=[];ci=[];return; end
%           if c > 2^(l-1), cl=[];cc=[];ci=[];return; end 
%           [dl,dc,di] = getChildren(obj,l,c);
%           [d1l,d1c,d1i] = getAllChildren(obj,dl(1),dc(1));
%           [d2l,d2c,d2i] = getAllChildren(obj,dl(2),dc(2));
%           cl = [dl d1l d2l];
%           cc = [dc d1c d2c];
%           ci = [di d1i d2i];
%        end