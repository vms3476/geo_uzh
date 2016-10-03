function [hdr] = readlas_hdr_function(fname)
% read out header information of las file -> Felix Morsdorf, UZH
  
hdr = mkhdr;  
fid = fopen(fname,'r','l');

hdr.FileSignature = char(fread(fid,4,'char')');
hdr.FileSourceID = fread(fid,1,'ushort');
hdr.Reserved = fread(fid,1,'ushort');
hdr.ProjectIDGUIDDATA1 = fread(fid,1,'ulong');
hdr.ProjectIDGUIDDATA2 = fread(fid,1,'ushort');
hdr.ProjectIDGUIDDATA3 = fread(fid,1,'ushort');
hdr.ProjectIDGUIDDATA4 = char(fread(fid,8,'char')');
hdr.VersionMajor = fread(fid,1,'char');
hdr.VersionMinor = fread(fid,1,'char');
hdr.SystemIdentifier = char(fread(fid,32,'char')');
hdr.GeneratingSoftware = char(fread(fid,32,'char')');
hdr.FileCreationDayofYear = fread(fid,1,'ushort');
hdr.FileCreationYear = fread(fid,1,'ushort');
hdr.HeaderSize = fread(fid,1,'ushort');
hdr.OffsetToPointData = fread(fid,1,'ulong');
hdr.NumberofVariableLengthRecords = fread(fid,1,'ulong');
hdr.PointDataFormatID = fread(fid,1,'char');
hdr.PointDataRecordLength = fread(fid,1,'ushort');
hdr.NumberOfPointRecords = fread(fid,1,'ulong');
hdr.NumberOfPointsbyReturn = fread(fid,5,'ulong');
hdr.XScaleFactor = fread(fid,1,'float64');
hdr.YScaleFactor = fread(fid,1,'float64');
hdr.ZScaleFactor = fread(fid,1,'float64');
hdr.XOffset = fread(fid,1,'float64');
hdr.YOffset = fread(fid,1,'float64');
hdr.ZOffset = fread(fid,1,'float64');
hdr.MaxX = fread(fid,1,'float64');
hdr.MinX = fread(fid,1,'float64');
hdr.MaxY = fread(fid,1,'float64');
hdr.MinY = fread(fid,1,'float64');
hdr.MaxZ = fread(fid,1,'float64');
hdr.MinZ = fread(fid,1,'float64');
fclose(fid); 
%----------------------------------------------------------------------
function hdr = mkhdr()
% set up blank hdr structure

hdr.FileSignature = [];
hdr.FileSourceID = [];
hdr.Reserved = [];
hdr.ProjectIDGUIDDATA1 = [];
hdr.ProjectIDGUIDDATA2 = [];
hdr.ProjectIDGUIDDATA3 = [];
hdr.ProjectIDGUIDDATA4 = [];
hdr.VersionMajor = [];
hdr.VersionMinor = [];
hdr.SystemIdentifier = [];
hdr.GeneratingSoftware = [];
hdr.FileCreationDayofYear = [];
hdr.FileCreationYear = [];
hdr.HeaderSize = [];
hdr.OffsetToPointData = [];
hdr.NumberofVariableLengthRecords = [];
hdr.PointDataFormatID = [];
hdr.PointDataRecordLength = [];
hdr.NumberOfPointRecords = [];
hdr.NumberOfPointsbyReturn = [];
hdr.XScaleFactor = [];
hdr.YScaleFactor = [];
hdr.ZScaleFactor = [];
hdr.XOffset = [];
hdr.YOffset = [];
hdr.ZOffset = [];
hdr.MaxX = [];
hdr.MinX = [];
hdr.MaxY = [];
hdr.MinY = [];
hdr.MaxZ = [];
hdr.MinZ = [];