function importWavefrontv2(obj,filename,path,mevis)
                 if nargin<4, mevis = false;end % Mevis = true, reading an .obj file as exported by MevisLab
                 if nargin<3, path = pwd; end
                 if ~strcmp(path(end),'/'), path = [path '/'];end
                 filename = [path filename];
                 if ~strcmp(filename(end-3:end),'.obj'), filename = [filename '.obj'];end
             % quick inreading of vertices, normals and texturecoordinates;
                 fid = fopen(filename,'r');
                 names = textscan(fid,'%s%*[^\n]');names = names{1};
                 fclose(fid);
                 fid = fopen(filename,'r');
                 data = textscan(fid,'%[^\n\r]');data = data{1};
                 fclose(fid);
             % reading in material lib file, and texturemaps
                 [UseMtl,Maps,startingpoints] = readMTL(filename);
             % reading Vertices
                 list = strcmp('v',names);index = find(list==1);
                 if isempty(index), return; end %nothing in the file
                 % small test to preallocate memory
                 Lyn = data{index(1)};tmp = sscanf(Lyn(2:end),'%f');
                 Location = nan*zeros(size(tmp,1),length(index));
                 for i=1:1:length(index)
                     Lyn = data{index(i)};Location(:,i) = sscanf(Lyn(2:end),'%f');
                 end
                 obj.Vertices = Location(1:3,:);
                 if size(Location,1)>3,  obj.VertexRGB = Location(4:end,:)./255; else, obj.VertexRGB = []; end 
                 clear list index Location;
             % reading normals 
                 % not required because computed in shape3D class on the fly    
             % reading texture coordinates
                 list = strcmp('vt',names);index = find(list==1);
                 if ~isempty(index)
                    uv = nan*zeros(2,length(index));
                    for i=1:1:length(index)
                        Lyn = data{index(i)};uv(:,i) = sscanf(Lyn(3:end),'%f');
                    end
                 else
                     uv = [];
                 end  
                 clear list index;
             % reading faces    
                list = strcmp('f',names);index = find(list==1);
                if ~isempty(index)
                    F3 = nan*zeros(3,length(index));F3t = F3;
                    F4 = nan*zeros(4,length(index));F4t = F4;
                    for i=1:1:length(index)
                        Lyn = data{index(i)};
                        Lyn=deblank(Lyn(3:end));
                    end
                end
end

function [UseMtl,Maps,startingpoints] = readMTL(filename)
        UseMtl = false;Maps = [];startingpoints = [];
        MTLName = [filename(1:end-4) '.mtl'];
        tmp = dir(MTLName);
        if isempty(tmp), return;end
        UseMtl = true;
        fid = fopen(MTLName,'r');
        MTL = textscan(fid,'%[^\n\r]');MTL = MTL{1};
        fclose(fid);
        nlines = length(MTL);
        % searching for mtl definition(s)
        newmtl = {};
        newmtlindex = [];
        for i=1:1:nlines
            str = MTL{i};
            if ~contains(str,'newmtl'), continue; end
            index = strfind(str,' ');
            newmtl{end+1} = str(index+1:end); %#ok<*AGROW>
            newmtlindex = [newmtlindex i];
        end
        nnewmtl = length(newmtl);
        % searching for map_Kd definitions
        map_Kd = cell(1,nnewmtl);
        for i=1:1:nnewmtl
            startindex = newmtlindex(i);
            if i<nnewmtl
               endindex = newmtlindex(i+1);
            else
               endindex = nlines;
            end
            for j=startindex:endindex
                str = MTL{j};
                if ~contains(str,'map_Kd'), continue; end
                index = strfind(str,' ');
                map_Kd{i} = str(index+1:end);
            end
        end
        nMaps = length(map_Kd);
        Maps = cell(1,nMaps);
        for i=1:1:nMaps
            if ~isempty(map_Kd{i}), Maps{i} = imread([path map_Kd{i}]);end
        end
        usemtl = +inf*ones(1,nnewmtl);
        list = strcmp('usemtl',names);index = find(list==1);
        listdata = data(index);
        for i=1:1:nnewmtl
            for j=1:1:length(listdata)
               if contains(listdata{j},newmtl{i})
                  usemtl(i) = index(j);
                  break;
               end
            end
        end
        
end

function [map,convertuv] = convert2SingleMap(Maps)




end
        
%                  if nargin<4, mevis = false;end % Mevis = true, reading an .obj file as exported by MevisLab
%                  if nargin<3, path = pwd; end
%                  if ~strcmp(path(end),'/'), path = [path '/'];end
%                  filename = [path filename];
%              % quick inreading of vertices, normals and texturecoordinates;
%                  fid = fopen(filename,'r');
%                  names = textscan(fid,'%s%*[^\n]');names = names{1};
%                  fclose(fid);
%                  fid = fopen(filename,'r');
%                  data = textscan(fid,'%[^\n\r]');
%                  data = data{1};
%                  fclose(fid);
%              % reading Vertices
%                  %disp('reading vertices');
%                  list = strcmp('v',names);index = find(list==1);
%                  if isempty(index), return; end %nothing in the file
%                  % small test to check whether rgb values are given
%                  Lyn = data{index(1)}; tmp = sscanf(Lyn(2:end),'%f');
%                  Location = nan*zeros(size(tmp,1),length(index));
%                  for i=1:1:length(index)
%                      Lyn = data{index(i)};
%                      Location(:,i) = sscanf(Lyn(2:end),'%f');
%                  end
%                  vertex = Location(1:3,:);
%                  if size(Location,1)>3 
%                      texturecolor = Location(4:end,:)./255; 
%                  else
%                      texturecolor = [];
%                  end
%                  clear list index Location;
%              % reading normals 
%                  % not required because computed in shape3D class on the fly    
%              % reading texture coordinates
%                  %disp('reading texture coordinates');
%                  list = strcmp('vt',names);index = find(list==1);
%                  if ~isempty(index)
%                     uv = nan*zeros(2,length(index));
%                     for i=1:1:length(index)
%                         Lyn = data{index(i)};
%                         uv(:,i) = sscanf(Lyn(3:end),'%f');
%                     end
%                  else
%                      uv = [];
%                  end  
%                  clear list index;
%              % reading faces
%                 list = strcmp('f',names);index = find(list==1);
%                 if ~isempty(index)
%                     F3 = nan*zeros(3,length(index));F3t = F3;
%                     F4 = nan*zeros(4,length(index));F4t = F4;
%                     for i=1:1:length(index)
%                         Lyn = data{index(i)};
%                         Lyn=deblank(Lyn(3:end));
%                         nvrts=length(findstr(Lyn,' ')); %#ok<*FSTR> % nr of vertices
%                         if ~mevis, nvrts = nvrts+1; end
%                         %nvrts=length(findstr(Lyn,' '))+1; %#ok<*FSTR> % nr of vertices
%                         if nvrts == 5, nvrts = 3; end
%                         if nvrts == 7, nvrts = 4; end
%                         nslash=length(findstr(Lyn,'/')); % nr of slashes
%                         switch nvrts
%                             case 3 % Triangles
%                                 switch nslash
%                                     case 0 % vertex
%                                         F3(:,i) = sscanf(Lyn,'%f');
%                                     case 3 % vertex\texture
%                                         f=sscanf(Lyn,'%f/%f');
%                                         F3t(:,i) = f([2 4 6]);
%                                         F3(:,i) = f([1 3 5]);
%                                     case 6 % vertex\texture\normal
%                                         f=sscanf(Lyn,'%f/%f/%f');
%                                         try
%                                          F3t(:,i) = f([2 5 8]);
%                                          F3(:,i) = f([1 4 7]);
%                                         catch % these are typically bad triangles, except for mevis files!!
%                                             if mevis
%                                                 f=sscanf(Lyn,'%f//%f');
%                                                 F3(:,i) = f([1 3 6]);
%                                             end
%                                         end
%                                     otherwise
%                                 end
%                             case 4 % Quadruples
%                                 switch nslash
%                                     case 0 % vertex
%         %                                 F4(:,i) = sscanf(Lyn,'%f');
%                                     case 4 % vertex\texture
%                                         f=sscanf(Lyn,'%f/%f');
%                                         F4t(:,i) = f([2 4 6 8]);
%                                         F4(:,i) = f([1 3 5 7]);
%                                     case 8 % vertex\texture\normal
%                                         f=sscanf(Lyn,'%f/%f/%f');
%                                         F4t(:,i) = f([2 5 8 11]);
%                                         F4(:,i) = f([1 4 7 10]);
%                                     otherwise
%                                 end
%                             otherwise
%                         end
%                     end
%                     % trimming the face lists, dropping the nan values
%                     index = find(~isnan(F3(1,:)));F3 = F3(:,index); %#ok<*FNDSB>
%                     index = find(~isnan(F3t(1,:)));F3t = F3t(:,index);
%                     index = find(~isnan(F4(1,:)));F4 = F4(:,index);
%                     index = find(~isnan(F4t(1,:)));F4t = F4t(:,index);
%                     %converting quads to triangles
%                     if ~isempty(F4)
%                         Tri1 = zeros(3,size(F4,2));
%                         Tri2 = zeros(3,size(F4,2));
%                                 for k=1:1:size(F4,2)
%                                     Tri1(:,k) = F4(1:3,k);
%                                     Tri2(1,k) = F4(1,k);
%                                     Tri2(2:end,k) = F4(3:end,k);
%                                 end
%                         F3 = [F3,Tri1,Tri2];
%                     end
%                     if ~isempty(F4t)
%                         Tri1 = zeros(3,size(F4t,2));
%                         Tri2 = zeros(3,size(F4t,2));
%                                 for k=1:1:size(F4t,2)
%                                     Tri1(:,k) = F4t(1:3,k);
%                                     Tri2(1,k) = F4t(1,k);
%                                     Tri2(2:end,k) = F4t(3:end,k);
%                                 end
%                         F3t = [F3t,Tri1,Tri2];
%                     end
%                 else
%                     F3 = [];F3t = F3;
%                     F4 = [];F4t = F4;
%                 end                
%                 
%              % Deleting non mesh points
%                 in_mesh = unique(F3(:))';
%                 not_nan = find(~isnan(vertex(1,:)));
%                 good = intersect(in_mesh,not_nan); 
%              % adjusting triangles
%                 [tmp,LOC] = ismember(F3,good);
%                 [~,j] = find(tmp==0);
%                 tmpindex = (1:size(F3,2));
%                 tmpindex = setdiff(tmpindex,j);
%                 F3 = LOC(:,tmpindex);
%              % adjusting texture coord information
%                 if ~isempty(F3t), F3t = F3t(:,tmpindex); end % this was commented before?
%              % adjusting vertex information    
%                 vertex = vertex(:,good);
%              %adjusting rgb information
%                 if ~isempty(texturecolor), texturecolor = texturecolor(:,good); end
%              % rearrangig texture coordinates
%                 if ~isempty(uv)           
%                    if ~isempty(F3t)
%                       [~,I,~] = unique(F3(:));
%                       T_index = F3t(I);
%                       uv = uv(:,T_index);    
%                    else
%                       uv = uv(:,good);
%                    end
%                 end
%               % Storing information in object
%                   obj.ColorMode = 'Single';
%                   obj.Vertices = vertex';
%                   obj.Faces = F3';
%                   if ~isempty(texturecolor)
%                       obj.VertexRGB = texturecolor';
%                       obj.ColorMode = 'Texture';
%                   end
%               % storing kind of wavefront in Userdatas
%               if ~mevis
%                   obj.UserData.ImportInfo = 'Wavefront Import';
%               else
%                   obj.UserData.ImportInfo = 'Mevislab Wavefront Import';
%               end
%               % reading texture file
%                   if isempty(uv), return; end % there is no texture image to look for
%                   name = filename(1:end-4);
%                   tfile = dir([name '.bmp']);
%                   if isempty(tfile),return; end% texture image does not exists
%                   uv(2,:) = 1-uv(2,:);
%                   obj.UV = uv;
%                   try
%                       im = imread([name '.bmp']);
%                       obj.TextureMap = im;
%                       obj.ColorMode = 'Texture';
%                   catch
%                   end
