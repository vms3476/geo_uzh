function [] = ers2hdr(filename,outfilename);
% function [] = ers2hdr(filename);
% Loads data header in ERMapper Format (*.ers)
% and converts to Envi HDR file   
%
% filename is XXXXX.ers - header file

%   
% Felix Morsdorf, RSL Zuerich, Jun 2005
  
  if nargin == 0
    [filename,pathname] = uigetfile('*.ers', ...
				    'Select Raw Data Header File ...');
    filename = [pathname,filename];
    k = findstr('.',filename)-1;
    outfilename = [filename(1:k)];    
  elseif nargin == 1
    k = findstr('.',filename)-1;
    outfilename = [filename(1:k)];
  end
  
  dfilename = [filename(1:k)]; 
  
  % load Header Information
  
  ers = readers(filename);
  
  % convert header information
  hdr = envi_header;
  
  hdr.description = ['{',ers.Version,' , ',ers.DataSetType,' converted by ers2hdr ', ...
                     datestr(clock),'}'];
  hdr.samples = num2str(ers.NrOfCellsPerLine);
  hdr.lines = num2str(ers.NrOfLines);
  hdr.bands = num2str(ers.NrOfBands);
  hdr.headeroffset = '0';
  hdr.filetype = 'Envi Standard';
  
  % parse data type
  if strcmp(upper(ers.CellType),upper('Unsigned8BitInteger'))
    hdr.datatype = '1';
  elseif strcmp(upper(ers.CellType),upper('Signed16BitInteger')) 
    hdr.datatype = '2';
  elseif strcmp(upper(ers.CellType),upper('Unsigned16BitInteger')) 
    hdr.datatype = '3';
  elseif strcmp(upper(ers.CellType),upper('IEEE4ByteReal'))
    hdr.datatype = '4';
  elseif strcmp(upper(ers.CellType),upper('IEEE8ByteReal'))
    hdr.datatype = '5';
  end
  hdr.interleave = 'bsq';
  hdr.sensortype = 'unknown';
  
  % parse byte swap
  if strcmp(upper(ers.ByteOrder),upper('MSBFirst'))
    hdr.byteorder = '1';
  elseif strcmp(upper(ers.ByteOrder),upper('LSBFirst'))
    hdr.byteorder = '0';
  end
  hdr.mapinfo = ['{',

  
  
  
  
%-------------------------------------------------------------------------------
function [hdr,vars] = envi_header();
  % initialize envi_header and set variables to defaults
  nvars = {'description','samples','lines','bands','header offset','file type', ...
	   'data type','interleave','sensor type','classes','class lookup', ...  
	   'class names','byte order','x start','y start','map info', ...  
	   'wavelength units','band names'};
  for i = 1:length(nvars)
    eval(['hdr.'n,deblank(nvars{i}),' = {}'])
  end
  
  
  
%-------------------------------------------------------------------------------  
function [success] = writehdr(hdr,file)
  % write envi header to file