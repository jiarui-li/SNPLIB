classdef morphableShape3D < superHandleClass
    properties
        Average;
        EigVec;
        EigVal;
        Neta = 1;
        Clamp = false;
        ClampValue = 2;
    end
    properties (Dependent = true)
        Dim;
        EigStd;
    end
    properties (Hidden = true)
        AlignmentTransformation;
    end
    properties (Hidden = true, Dependent = true)
        AverageVEC;
    end
    methods %CONSTRUCTOR
        function obj = morphableShape3D(varargin)
            obj = obj@superHandleClass(varargin{:});
        end
    end 
    methods %GETTING
        function out = get.AverageVEC(obj)
            if isempty(obj.Average), out = []; return; end
            out = obj.Average.Vertices';out = out(:);
        end
        function out = get.Average(obj)
            out = obj.Average;
            if ~superHandleClass.isH(out), out = []; end
        end
        function out = get.Dim(obj)
            if isempty(obj.EigVal), out = []; return;end
            out = length(obj.EigVal);
        end
        function out = get.EigStd(obj)
            out = sqrt(obj.EigVal);
        end
    end
    methods %SETTING
    end
    methods %INTERFACING
        function out = alignIntoSpace(obj,in,w)
            if nargin<3, w = ones(1,in.nVertices);end
            obj.AlignmentTransformation = morphableShape3D.getTransformation(obj.Average,in,true,w); 
            out = morphableShape3D.evalTransformation(obj.AlignmentTransformation,in);  
        end
        function out = projectIntoSpace(obj,in,w)
            if nargin<3, w = ones(1,in.nVertices);end
            % converting input to Vector representation
            inVEC = in.Vertices';
            inVEC = double(inVEC(:));
            inVEC = inVEC-obj.AverageVEC;
            wVEC = repmat(w,3,1);
            wVEC = double(wVEC(:));
            % setting the regularisation
            sigma = obj.Neta*obj.EigVal(1);
            % getting weight matrix and perform projection
            A = spdiags(wVEC,0,speye(length(inVEC),length(inVEC)));
            AQ = A*obj.EigVec*diag(obj.EigVal);
            [U,W,V] = svd(AQ,'econ');
            W = diag(W);
            W = W./((W.^2)+ones(size(W))*sigma);
            out = diag(obj.EigVal)*V*diag(W)*U'*A*inVEC;
        end
        function out = getShapeInstance(obj,scores)
           if isempty(obj.Average), out = []; return; end
           VEC = obj.AverageVEC + obj.EigVec*scores;
           out = clone(obj.Average);
           out.Vertices = reshape(VEC,3,obj.Average.nVertices)';
        end
        function [out,scores] = compute_morphable_transformation(obj,in,w)
            if nargin<3, w = ones(1,in.nVertices);end
            % step 1, align with average
            out = alignIntoSpace(obj,in,w);
            % step 2, project into PCA space
            out = projectIntoSpace(obj,out,w);
            scores = out;
            % step 3, clamping
            if obj.Clamp, out = sign(out).*min([abs(out) abs(obj.EigStd).*obj.ClampValue],[],2); end
            % step 4, recreate a shape
            out = getShapeInstance(obj,out);
            % step 5, transform back to original shape coordinate system
            out = morphableShape3D.evalInverseTransformation(obj.AlignmentTransformation,out);
            out.Vertices = double(out.Vertices);
        end   
    end
    methods (Static = true)
        function T = getTransformation(Target,Floating,scale,w)
             q = Target.Vertices';
             p = Floating.Vertices';
             nbpts = size(p,2);
             index = find(w);% find points not having weights == 0;
             if length(index)<nbpts
                 p = p(:,index);
                 q = q(:,index);
                 w = w(:,index);
             end
             totalw = sum(w,2);
             rw = repmat(w,3,1);
             nbpts = size(p,2);
             % align centers
             centerq = sum(rw.*q,2)./repmat(totalw,3,1);
             centerp = sum(rw.*p,2)./repmat(totalw,3,1);
             q = q-repmat(centerq,1,nbpts);
             p = p-repmat(centerp,1,nbpts);
             % Scale both by making mean length 1
             if scale
                 lengths = sqrt(sum((rw.*p).^2,1));
                 meanLength1 = sum(lengths,2)/totalw;
                 p = p./meanLength1;
                 lengths = sqrt(sum((rw.*q).^2,1));
                 meanLength2 = sum(lengths,2)/totalw;
                 q = q./meanLength2;
                 scalefactor = meanLength2/meanLength1;
             else
                 scalefactor = 1;% keep scale fixed
             end
             % Rotate
             [U,S,V] = svd(rw.*q*p');
             H = V*sign(S)*U';
             H = H';
             % putting it all together
             T.Scale=scalefactor;
             T.Rotation = H;
             transform = T.Scale*H;
             Ta=eye(4);
             Ta(1:3,4)=centerq;
             Tb=eye(4);
             Tb(1:3,4)=-centerp;
             R=eye(4);
             R(1:3,1:3)=transform;
             Tout=Ta*R*Tb;
             T.Translation = Tout(1:3,4)/T.Scale;
        end
        function out = evalTransformation(T,in)
             out = clone(in);
             out.Vertices = (T.Scale*(T.Rotation*(out.Vertices')+repmat(T.Translation,1,out.nVertices)))';
        end
        function out = evalInverseTransformation(T,in)
             out = clone(in);
             invRotation = T.Rotation^-1;
             out.Vertices = (invRotation*(out.Vertices'/T.Scale-repmat(T.Translation,1,out.nVertices)))';
        end
    end
end