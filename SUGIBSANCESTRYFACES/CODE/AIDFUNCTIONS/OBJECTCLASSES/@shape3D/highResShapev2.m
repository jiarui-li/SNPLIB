function out = highResShapev2(obj,TemplateTextureInfo)
% the idea is to create another high resolution shape instance, equal to
% the resolution of the TextureMap, to obtain a photorealistic rendering in
% matlab

    out = shape3D;% create empty shape
    tri3D = triangulation(TemplateTextureInfo.tri3D.ConnectivityList,obj.Vertices);% Triangulation specific for obj shape
    out.Vertices = barycentricToCartesian(tri3D,double(TemplateTextureInfo.FIndex),TemplateTextureInfo.BAR);% compute vertices from pixels
    
    triUV = triangulation(TemplateTextureInfo.tri3D.ConnectivityList,[obj.UV zeros(obj.nVertices,1)]);% Triangulation specific for obj shape
    out.UV = barycentricToCartesian(triUV,double(TemplateTextureInfo.FIndex),TemplateTextureInfo.BAR);% compute vertices from pixels
    out.UV = out.UV(:,1:2);
    
    out.TextureMap = obj.TextureMap;
    
    out.Faces = TemplateTextureInfo.ImageTri.ConnectivityList;% assign precomputation triangulation on all pixels
    
%     R = squeeze(obj.TextureMap(:,:,1));% extract red channel
%     G = squeeze(obj.TextureMap(:,:,2));% extract green channel
%     B = squeeze(obj.TextureMap(:,:,3));% extract blue channel
%     out.VertexRGB = [R(:) G(:) B(:)];% assign colors to shape object

end