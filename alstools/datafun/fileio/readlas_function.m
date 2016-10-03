function [data] = readlas_function(fname)
% read out header & data las file -> Felix Morsdorf, UZH, 2011 
hdr = readlas_hdr_function(fname);
[fid] = fopen(fname,'r','l');
[data] = rekpoints(fid,hdr);
fclose(fid);
end

%-------------------------------------------------------------------------
function [raw] = rekpoints(fid,hdr)
% function to reconstruct point data from point record   
num = hdr.NumberOfPointRecords;
x = ones(1,num);raw.y = x;raw.z = x;raw.int = x;
raw.x = x;
  
if 1
	raw.rnnr = x;raw.nrrt =x;
    raw.ReturnNumber = x; raw.NumberOfReturns = x;
    raw.ScanDirection = x;raw.EdgeFlightLine = x;
    raw.Classification = x; raw.ScanAngle = x;
    raw.UserData = x; raw.PointSourceID = x;
    raw.GPSTime = x;
end
  
clear x;
  
skip = hdr.PointDataRecordLength;
fseek(fid,hdr.OffsetToPointData,-1);
  
raw.x = fread(fid,num,'long',skip-4);
fseek(fid,hdr.OffsetToPointData,-1);fseek(fid,4,0);
  
raw.y = fread(fid,num,'long',skip-4);
fseek(fid,hdr.OffsetToPointData,-1);fseek(fid,8,0);
  
raw.z = fread(fid,num,'long',skip-4);
fseek(fid,hdr.OffsetToPointData,-1);fseek(fid,12,0);
  
raw.int = fread(fid,num,'ushort',skip-2);
  
if 1
	fseek(fid,hdr.OffsetToPointData,-1);fseek(fid,14,0);
    % return number and number of returns are coded with 3 bits each!
    raw.rnnr = fread(fid,num,'*ubit6',((skip-1)*8)+2);

    fseek(fid,hdr.OffsetToPointData,-1);fseek(fid,15,0);
    raw.Classification = fread(fid,num,'char',skip-1);
    
    fseek(fid,hdr.OffsetToPointData,-1);fseek(fid,16,0);
    raw.ScanAngle = fread(fid,num,'char=>uchar',skip-1);
    
    fseek(fid,hdr.OffsetToPointData,-1);fseek(fid,17,0);
    raw.UserData = fread(fid,num,'char',skip-1);
    
    fseek(fid,hdr.OffsetToPointData,-1);fseek(fid,18,0);
    raw.PointSourceID = fread(fid,num,'ushort',skip-2);
    
    fseek(fid,hdr.OffsetToPointData,-1);fseek(fid,20,0);
    raw.GPSTime = fread(fid,num,'double',skip-8);
end
  
% correct for missing GPSTime values
diffl = length(raw.GPSTime) - length(raw.x);
if diffl ~= 0
	sampl = mean(diff(raw.GPSTime));
    for i = 1:diffl
    	raw.GPSTime(end+i) = raw.GPSTime(end+i-1)+sampl;
    end
end
  
raw.x = (raw.x * hdr.XScaleFactor) + hdr.XOffset;
raw.y = (raw.y * hdr.YScaleFactor) + hdr.YOffset;
raw.z = (raw.z * hdr.ZScaleFactor) + hdr.ZOffset;
  
end