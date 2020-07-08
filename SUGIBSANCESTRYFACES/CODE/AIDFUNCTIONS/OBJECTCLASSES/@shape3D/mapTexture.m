function out = mapTexture(obj,Template,TemplateTextureInfo,origmesh)
% specific function, not aimed to be distributed
% To aid in texture extraction, this function prepares information on a
% template (obj needs to be a template or refscan as such)

% The idea is as follows: For each pixel in the mapped texturemap, go to 3D
% using the mapped shape triangulation, then find the closest three points
% on the original shape and retrieve UV coordinates in teh original shape
% Texture map, by interpolating the UV coordinates from the three closest
% points.
warning off;
    obj.TextureMap = Template.TextureMap;% Preassign the Template's texturemap to the object
    obj.UV = Template.UV;% Transfer the templates UV coordinates, as these are standard to use.
 
    tri3D = triangulation(TemplateTextureInfo.tri3D.ConnectivityList,obj.Vertices);% define objects triangulation for points defined on the 2D texturemap
    % location of each pixel onto the given face in obj and origobj.
    CART = barycentricToCartesian(tri3D,double(TemplateTextureInfo.FIndex),TemplateTextureInfo.BAR);% For each pixel compute its position in 3D
    % getting the bary coordinates on the original mesh
    [bar,index] = cart2baryKNN(origmesh.Vertices,CART);% For each pixel in 3D, what are the three closest points on the original mesh and express in barycentric coordinates
    clear CART;% clear some memory
    n = size(bar,1);
    UV = zeros(n,2);
    origuv = origmesh.UV;
    parfor i=1:1:n% for each pixel in 3D, find its UV position on the originalMesh texturemap 2D domain.
        warning off;
        %i=100;
        a = origuv(index(i,1),:); %#ok<*PFBNS>
        b = origuv(index(i,2),:);
        c = origuv(index(i,3),:);
        coherenceWeight = zeros(3,1);
        if sqrt(sum((a-b).^2)) < 0.01
           coherenceWeight(1)=1;
           coherenceWeight(2)=1;
        end
        if sqrt(sum((a-c).^2)) < 0.01
           coherenceWeight(1)=1;
           coherenceWeight(3)=1;
        end
        if sqrt(sum((c-b).^2)) < 0.01
           coherenceWeight(3)=1;
           coherenceWeight(2)=1;
        end
        w1 = coherenceWeight(1).*bar(i,1);
        w2 = coherenceWeight(2).*bar(i,2);
        w3 = coherenceWeight(3).*bar(i,3);
        sumw = w1+w2+w3;w1 = w1/sumw;w2 = w2/sumw;w3 = w3/sumw;
        if sumw>0.001 % at least to coherent UV coordinates
           tmpuv = w1.*a+w2.*b+w3.*c;
        else % select UV of closest point
           tmpuv = a;
        end 
        intmp = inpolygon(tmpuv(1),tmpuv(2),[a(1);b(1);c(1)],[a(2);b(2);c(2)]);      
        if ~intmp, tmpuv = a;end% select UV closest point, when projected point is not within the triangle of interpolating points (indicates bad interpolation)
        tmpuv(tmpuv<0) = 0;
        tmpuv(tmpuv>1) = 1;
        UV(i,:) = tmpuv;
    end
    out = textureMapValues(origmesh,UV);% for each pixel, retrieve the RGB values from the originalMesh textureMap, using the established UV coordinates
    obj.TextureMap = reshape(out,size(obj.TextureMap));% assign the RGB values for each pixel into an Image format
warning on;    
%     figure;scatter(UV(:,1),UV(:,2),10,inok,'filled');    
%    
%     figure;plot(UV(:,1),UV(:,2),'b.');
%     
%     figure;plot(origuv(:,1),-1*origuv(:,2),'r.');
%     
%     
%     figure;imshow(obj.TextureMap)
%     
%     cartShape = shape3D;
%     cartShape.Vertices = CART;
%     cartShape.VertexRGB = out;
%     v = viewer(cartShape);
% %     viewer(obj,v);
% %     

end

