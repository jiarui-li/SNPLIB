function exportWavefrontv1(obj,filename,path)
    % changing directory
         if nargin<3, path = pwd; end
         if ~strcmp(path(end),'/'), path = [path '/'];end
         filename = [path filename];
         if strcmp(filename(end-3:end),'.obj'),filename = filename(1:end-4);end
     % Opening file and write Heading
         str = [filename '.obj'];
         fid = fopen(str,'w');
         fprintf(fid,['# OBJ exported by ' obj.Type ' matlab class \n']);
         fprintf(fid,'\n');
         str = ['mtllib ' filename '.mtl'];
         fprintf(fid,str);
         fprintf(fid,'\n');
         fprintf(fid,'\n');
         fprintf(fid, '# %d vertex\n', obj.nVertices);
         fprintf(fid,'\n');
     % Write according to Texture information
         if isempty(obj.TextureMap)||isempty(obj.UV)% No texture Map
             % writing the vertex information            
             if ~isempty(obj.VertexRGB)% RGB Colours per Vertex
                fprintf(fid,'v %f %f %f %f %f %f\n',[obj.Vertices';obj.VertexRGB*255']);    
             else
                fprintf(fid,'v %f %f %f\n',obj.Vertices'); 
             end
             fprintf(fid,'\n');
             str = ['g ' filename];fprintf(fid,str);
             fprintf(fid,'\n');
             str = ['usemtl ' filename];fprintf(fid,str);
             fprintf(fid,'\n');
             str = ['s'];fprintf(fid,str);
             fprintf(fid,'\n');
             fprintf(fid,'\n');
             % writing face information
             fprintf(fid, '# %d faces\n', obj.nFaces);
             fprintf(fid,'\n');
             fprintf(fid,'f %.0f %.0f %.0f\n',obj.Faces');
             fclose(fid);
         else % With TextureMap
             % writing image
             warning off;
             imwrite(obj.TextureMap,[filename '.bmp'],'bmp');
             warning on;
             UV = obj.UV';
             UV(2,:) = 1-UV(2,:);
             % writing vertices
             fprintf(fid,'v %f %f %f\n',obj.Vertices');
             fprintf(fid,'\n');
             % writing Texture Coordinates
             fprintf(fid,'vt %f %f\n',UV);
             clear UV;
             fprintf(fid,'\n');
             str = ['g ' filename];fprintf(fid,str);
             fprintf(fid,'\n');
             str = ['usemtl ' filename];fprintf(fid,str);
             fprintf(fid,'\n');
             str = ['s'];fprintf(fid,str); %#ok<*NBRAK>
             fprintf(fid,'\n');
             fprintf(fid,'\n');
             Tri = obj.Faces';
             tmp = zeros(6,size(Tri,2));
             tmp(1,:) = Tri(1,:);tmp(2,:) = Tri(1,:);
             tmp(3,:) = Tri(2,:);tmp(4,:) = Tri(2,:);
             tmp(5,:) = Tri(3,:);tmp(6,:) = Tri(3,:);
             fprintf(fid,'f %.0f/%.0f %.0f/%.0f %.0f/%.0f\n',tmp);
             fclose(fid);
         end
     % Saving the relectance information (.mtl file)
        str = [filename '.mtl'];
        fid = fopen(str,'w');
        str = ['newmtl ' filename];fprintf(fid,str);fprintf(fid,'\n');
        str = ['ka 0.3 0.3 0.3'];fprintf(fid,str);fprintf(fid,'\n');
        str = ['kd 0.8 0.8 0.8'];fprintf(fid,str);fprintf(fid,'\n');
        str = ['ks 0 0 0'];fprintf(fid,str);fprintf(fid,'\n');
        str = ['illum 0'];fprintf(fid,str);fprintf(fid,'\n');
        str = ['map_Ka ' filename '.bmp'];fprintf(fid,str);fprintf(fid,'\n');
        str = ['map_Kd ' filename '.bmp'];fprintf(fid,str);fprintf(fid,'\n');
        str = ['map_Ks ' filename '.bmp'];fprintf(fid,str);fprintf(fid,'\n');
        fclose(fid);   
end