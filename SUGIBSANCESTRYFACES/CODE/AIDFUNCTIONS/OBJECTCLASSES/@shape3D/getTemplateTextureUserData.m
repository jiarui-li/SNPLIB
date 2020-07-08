function out = getTemplateTextureUserData(obj,res)
% specific function, not aimed to be distributed
% To aid in texture extraction, this function prepares information on a
% template (obj needs to be a template or refscan as such)

    randIm = randi(255,res,res,3,'uint8');
    obj.TextureMap = randIm;
    [X,Y] = textureMapGrid(obj);   
    QP = [X(:),Y(:)];
    nQP = size(QP,1);
    
    out.ImageTri = delaunayTriangulation(QP);
    tri2D = delaunayTriangulation(obj.UV);out.tri2D = tri2D;
    tri3D = triangulation(tri2D.ConnectivityList, obj.Vertices);out.tri3D = tri3D;
    
    F = getFaceIndexTri2D(QP,tri2D);
    out.FIndex = superHandleClass.convertUInt(F);
    
    BAR = cartesianToBarycentric(tri2D,F,QP);out.BAR = BAR;
    CART = barycentricToCartesian(tri3D,F,BAR);out.CART = CART;

    cartShape = shape3D;
    cartShape.Vertices = CART;
    v = viewer(cartShape);
    viewer(obj,v);
    
end

function F = getFaceIndexTri2D(QP,tri)
         nQP = size(QP,1);
         nF = size(tri.ConnectivityList,1);
         F = zeros(nQP,1);  
         for f=1:nF
             face = tri.ConnectivityList(f,:);
             IN = inpolygon(QP(:,1),QP(:,2),tri.Points(face,1),tri.Points(face,2));
             F(IN) = f;
         end
end