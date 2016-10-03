function[]=mat2dxf(fname,vertex,faces) 
%
% Given a filename and a the vertex and faces of the 3D model that you want
% to convert to DXF
% for more information email: jonas_van_den_bergh@hotmail.com
% try at your own risk, refer also to writedxf.m from Greg Siegle 
%
name=sprintf('%s.dxf',fname);
fid=fopen(name,'w');
fprintf(fid,'999\ncreated by Matlab\n0\nSECTION\n2\nENTITIES\n0\n');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[aantalVertex,zever] = size(vertex(:,1));
[aantalFaces,zever] = size(faces(:,1));

fprintf(fid,'POLYLINE\n8\n%s\n66\n1\n70\n64\n71\n%d\n72\n%d\n62\n144\n0\n',fname,aantalVertex, (2*aantalFaces));

%vertex
for a=1:aantalVertex
      %header
      fprintf(fid,'VERTEX\n8\n%s\n',fname);
      
      %coordinates
      fprintf(fid,'10\n%.4f\n20\n%.4f\n30\n%.4f\n',vertex(a,1),vertex(a,2),vertex(a,3));
      
      %end
      fprintf(fid,'70\n192\n0\n');      
      
end
if 0
%face
for a=1:aantalFaces
      %header
      fprintf(fid,'VERTEX\n8\n%s\n10\n0\n20\n0\n30\n0\n70\n128\n',fname);
      
      %connections
      fprintf(fid,'71\n%d\n72\n%d\n73\n-%d\n',faces(a,1),faces(a,2),faces(a,3));
      
      %end
      fprintf(fid,'0\n');
      
      %header
      fprintf(fid,'VERTEX\n8\n%s\n10\n0\n20\n0\n30\n0\n70\n128\n',fname);
      
      %connections
      fprintf(fid,'71\n%d\n72\n%d\n73\n-%d\n',faces(a,3),faces(a,4),faces(a,1));
            
      %end
      fprintf(fid,'0\n');      
      
end
end
%close 
fprintf(fid,'SEQEND\n8\n%s\n0\n',fname);
fprintf(fid,'ENDSEC\n0\nEOF\n');
fclose(fid);