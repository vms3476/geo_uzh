function geotiffwrite_function(filename, bbox, image, bit_depth, option)
%                                                    (Version 29august2011)
%--------------------------------------------------------------------------
   % set default inputs
   if nargin < 5
      option = [];
      if isempty(bbox)
         help geotiffwrite;
         error('Because you did not specify bounding box, please assign proper values to option argument');
      end
   end

   if nargin < 4
      bit_depth = [];
   end

   if nargin < 3
      error('Usage: geotiffwrite(filename, bbox, image, [bit_depth, option])');
   end

   if isempty(bit_depth)
      bit_depth = 32;
   end

   if isfield(option, 'ColorMap') && ( abs(bit_depth)==32 || abs(bit_depth)==1 )
      error('If you set values to option.ColorMap, you cannot write in 32-bit depth or 1-bit depth.');
   end

   if isfield(option, 'NaN')
      image(isnan(image)) = option.NaN;
   end

   % set Orientation, 274
   if isfield(option, 'Orientation')
      ifd.Orientation = option.Orientation;

      for i = 1:size(image,3)
         switch ifd.Orientation
         case 2
            image(:,:,i) = fliplr(image(:,:,i));
         case 3
            image(:,:,i) = fliplr(flipud(image(:,:,i)));
         case 4
            image(:,:,i) = flipud(image(:,:,i));
         case 5
            image(:,:,i) = flipud(rot90(image(:,:,i)));
         case 6
            image(:,:,i) = rot90(image(:,:,i));
         case 7
            image(:,:,i) = fliplr(rot90(image(:,:,i)));
         case 8
            image(:,:,i) = fliplr(flipud(rot90(image(:,:,i))));
         end
      end
   else
      ifd.Orientation = 1;
   end

   if abs(bit_depth) == 1 && mod(size(image,2),8) ~= 0
      error('Please use integral multiple of 8 for 1-bit image size.');
   end

   GeoAsciiParamsTag = '';
   GeoDoubleParamsTag = [];
   GeoAsciiOffset = 0;
   GeoDoubleOffset = 0;
   NumberOfKeys = 0;
   GeoKeyDirectoryTag = [1 1 0 NumberOfKeys];

   % set GTModelTypeGeoKey, 1024
   NumberOfKeys = NumberOfKeys + 1;
   GeoKeyDirectoryTag(1, 4) = NumberOfKeys;

   if isfield(option, 'GTModelTypeGeoKey') && option.GTModelTypeGeoKey ~= 2
      code = option.GTModelTypeGeoKey;
      bbox = [];
   else
      code = 2;					% ModelTypeGeographic
   end

   GeoKeyDirectoryTag = [GeoKeyDirectoryTag; [1024 0 1 code]];

   % set GTRasterTypeGeoKey, 1025
   NumberOfKeys = NumberOfKeys + 1;
   GeoKeyDirectoryTag(1, 4) = NumberOfKeys;

   if isfield(option, 'GTRasterTypeGeoKey')
      code = option.GTRasterTypeGeoKey;
   else
      code = 1;					% RasterPixelIsArea
   end

   gifd.GTRasterTypeGeoKey = code;
   GeoKeyDirectoryTag = [GeoKeyDirectoryTag; [1025 0 1 code]];

   %  set GTCitationGeoKey, 1026
   if isfield(option, 'GTCitationGeoKey')
      NumberOfKeys = NumberOfKeys + 1;
      GeoKeyDirectoryTag(1, 4) = NumberOfKeys;
      code = [option.GTCitationGeoKey(:); '|']; cnt = length(code);
      GeoKeyDirectoryTag = [GeoKeyDirectoryTag; [1026 34737 cnt GeoAsciiOffset]];
      GeoAsciiParamsTag = [GeoAsciiParamsTag; code];
      GeoAsciiOffset = GeoAsciiOffset + cnt;
   end

   % set GeographicTypeGeoKey, 2048
   if isfield(option, 'GeographicTypeGeoKey')
      NumberOfKeys = NumberOfKeys + 1;
      GeoKeyDirectoryTag(1, 4) = NumberOfKeys;
      code = option.GeographicTypeGeoKey;
      GeoKeyDirectoryTag = [GeoKeyDirectoryTag; [2048 0 1 code]];
   elseif ~isempty(bbox)
      NumberOfKeys = NumberOfKeys + 1;
      GeoKeyDirectoryTag(1, 4) = NumberOfKeys;
      code = 4326;				% GCS_WGS_84
      GeoKeyDirectoryTag = [GeoKeyDirectoryTag; [2048 0 1 code]];
   end

   % set GeogCitationGeoKey, 2049
   if isfield(option, 'GeogCitationGeoKey')
      NumberOfKeys = NumberOfKeys + 1;
      GeoKeyDirectoryTag(1, 4) = NumberOfKeys;
      code = [option.GeogCitationGeoKey(:); '|']; cnt = length(code);
      GeoKeyDirectoryTag = [GeoKeyDirectoryTag; [2049 34737 cnt GeoAsciiOffset]];
      GeoAsciiParamsTag = [GeoAsciiParamsTag; code];
      GeoAsciiOffset = GeoAsciiOffset + cnt;
   end

   % set GeogGeodeticDatumGeoKey, 2050
   if isfield(option, 'GeogGeodeticDatumGeoKey')
      NumberOfKeys = NumberOfKeys + 1;
      GeoKeyDirectoryTag(1, 4) = NumberOfKeys;
      code = option.GeogGeodeticDatumGeoKey;
      GeoKeyDirectoryTag = [GeoKeyDirectoryTag; [2050 0 1 code]];
   end

   % set GeogPrimeMeridianGeoKey, 2051
   if isfield(option, 'GeogPrimeMeridianGeoKey')
      NumberOfKeys = NumberOfKeys + 1;
      GeoKeyDirectoryTag(1, 4) = NumberOfKeys;
      code = option.GeogPrimeMeridianGeoKey;
      GeoKeyDirectoryTag = [GeoKeyDirectoryTag; [2051 0 1 code]];
   end

   % set GeogLinearUnitsGeoKey, 2052
   if isfield(option, 'GeogLinearUnitsGeoKey')
      NumberOfKeys = NumberOfKeys + 1;
      GeoKeyDirectoryTag(1, 4) = NumberOfKeys;
      code = option.GeogLinearUnitsGeoKey;
      GeoKeyDirectoryTag = [GeoKeyDirectoryTag; [2052 0 1 code]];
   end

   % set GeogLinearUnitSizeGeoKey, 2053
   if isfield(option, 'GeogLinearUnitSizeGeoKey')
      NumberOfKeys = NumberOfKeys + 1;
      GeoKeyDirectoryTag(1, 4) = NumberOfKeys;
      code = option.GeogLinearUnitSizeGeoKey;
      GeoKeyDirectoryTag = [GeoKeyDirectoryTag; [2053 34736 1 GeoDoubleOffset]];
      GeoDoubleParamsTag = [GeoDoubleParamsTag; code];
      GeoDoubleOffset = GeoDoubleOffset + 1;
   end

   % set GeogAngularUnitsGeoKey, 2054
   if isfield(option, 'GeogAngularUnitsGeoKey')
      NumberOfKeys = NumberOfKeys + 1;
      GeoKeyDirectoryTag(1, 4) = NumberOfKeys;
      code = option.GeogAngularUnitsGeoKey;
      GeoKeyDirectoryTag = [GeoKeyDirectoryTag; [2054 0 1 code]];
   end

   % set GeogAngularUnitSizeGeoKey, 2055
   if isfield(option, 'GeogAngularUnitSizeGeoKey')
      NumberOfKeys = NumberOfKeys + 1;
      GeoKeyDirectoryTag(1, 4) = NumberOfKeys;
      code = option.GeogAngularUnitSizeGeoKey;
      GeoKeyDirectoryTag = [GeoKeyDirectoryTag; [2055 34736 1 GeoDoubleOffset]];
      GeoDoubleParamsTag = [GeoDoubleParamsTag; code];
      GeoDoubleOffset = GeoDoubleOffset + 1;
   end

   % set GeogEllipsoidGeoKey, 2056
   if isfield(option, 'GeogEllipsoidGeoKey')
      NumberOfKeys = NumberOfKeys + 1;
      GeoKeyDirectoryTag(1, 4) = NumberOfKeys;
      code = option.GeogEllipsoidGeoKey;
      GeoKeyDirectoryTag = [GeoKeyDirectoryTag; [2056 0 1 code]];
   end

   % set GeogSemiMajorAxisGeoKey, 2057
   if isfield(option, 'GeogSemiMajorAxisGeoKey')
      NumberOfKeys = NumberOfKeys + 1;
      GeoKeyDirectoryTag(1, 4) = NumberOfKeys;
      code = option.GeogSemiMajorAxisGeoKey;
      GeoKeyDirectoryTag = [GeoKeyDirectoryTag; [2057 34736 1 GeoDoubleOffset]];
      GeoDoubleParamsTag = [GeoDoubleParamsTag; code];
      GeoDoubleOffset = GeoDoubleOffset + 1;
   end

   % set GeogSemiMinorAxisGeoKey, 2058
   if isfield(option, 'GeogSemiMinorAxisGeoKey')
      NumberOfKeys = NumberOfKeys + 1;
      GeoKeyDirectoryTag(1, 4) = NumberOfKeys;
      code = option.GeogSemiMinorAxisGeoKey;
      GeoKeyDirectoryTag = [GeoKeyDirectoryTag; [2058 34736 1 GeoDoubleOffset]];
      GeoDoubleParamsTag = [GeoDoubleParamsTag; code];
      GeoDoubleOffset = GeoDoubleOffset + 1;
   end

   % set GeogInvFlatteningGeoKey, 2059
   if isfield(option, 'GeogInvFlatteningGeoKey')
      NumberOfKeys = NumberOfKeys + 1;
      GeoKeyDirectoryTag(1, 4) = NumberOfKeys;
      code = option.GeogInvFlatteningGeoKey;
      GeoKeyDirectoryTag = [GeoKeyDirectoryTag; [2059 34736 1 GeoDoubleOffset]];
      GeoDoubleParamsTag = [GeoDoubleParamsTag; code];
      GeoDoubleOffset = GeoDoubleOffset + 1;
   end

   % set GeogAzimuthUnitsGeoKey, 2060
   if isfield(option, 'GeogAzimuthUnitsGeoKey')
      NumberOfKeys = NumberOfKeys + 1;
      GeoKeyDirectoryTag(1, 4) = NumberOfKeys;
      code = option.GeogAzimuthUnitsGeoKey;
      GeoKeyDirectoryTag = [GeoKeyDirectoryTag; [2060 0 1 code]];
   end

   % set GeogPrimeMeridianLongGeoKey, 2061
   if isfield(option, 'GeogPrimeMeridianLongGeoKey')
      NumberOfKeys = NumberOfKeys + 1;
      GeoKeyDirectoryTag(1, 4) = NumberOfKeys;
      code = option.GeogPrimeMeridianLongGeoKey;
      GeoKeyDirectoryTag = [GeoKeyDirectoryTag; [2061 34736 1 GeoDoubleOffset]];
      GeoDoubleParamsTag = [GeoDoubleParamsTag; code];
      GeoDoubleOffset = GeoDoubleOffset + 1;
   end

   % set ProjectedCSTypeGeoKey, 3072
   if isfield(option, 'ProjectedCSTypeGeoKey')
      NumberOfKeys = NumberOfKeys + 1;
      GeoKeyDirectoryTag(1, 4) = NumberOfKeys;
      code = option.ProjectedCSTypeGeoKey;
      GeoKeyDirectoryTag = [GeoKeyDirectoryTag; [3072 0 1 code]];
   end

   % set PCSCitationGeoKey, 3073
   if isfield(option, 'PCSCitationGeoKey')
      NumberOfKeys = NumberOfKeys + 1;
      GeoKeyDirectoryTag(1, 4) = NumberOfKeys;
      code = [option.PCSCitationGeoKey(:); '|']; cnt = length(code);
      GeoKeyDirectoryTag = [GeoKeyDirectoryTag; [3073 34737 cnt GeoAsciiOffset]];
      GeoAsciiParamsTag = [GeoAsciiParamsTag; code];
      GeoAsciiOffset = GeoAsciiOffset + cnt;
   end

   % set ProjectionGeoKey, 3074
   if isfield(option, 'ProjectionGeoKey')
      NumberOfKeys = NumberOfKeys + 1;
      GeoKeyDirectoryTag(1, 4) = NumberOfKeys;
      code = option.ProjectionGeoKey;
      GeoKeyDirectoryTag = [GeoKeyDirectoryTag; [3074 0 1 code]];
   end

   % set ProjCoordTransGeoKey, 3075
   if isfield(option, 'ProjCoordTransGeoKey')
      NumberOfKeys = NumberOfKeys + 1;
      GeoKeyDirectoryTag(1, 4) = NumberOfKeys;
      code = option.ProjCoordTransGeoKey;
      GeoKeyDirectoryTag = [GeoKeyDirectoryTag; [3075 0 1 code]];
   end

   % set ProjLinearUnitsGeoKey, 3076
   if isfield(option, 'ProjLinearUnitsGeoKey')
      NumberOfKeys = NumberOfKeys + 1;
      GeoKeyDirectoryTag(1, 4) = NumberOfKeys;
      code = option.ProjLinearUnitsGeoKey;
      GeoKeyDirectoryTag = [GeoKeyDirectoryTag; [3076 0 1 code]];
   end

   % set ProjLinearUnitSizeGeoKey, 3077
   if isfield(option, 'ProjLinearUnitSizeGeoKey')
      NumberOfKeys = NumberOfKeys + 1;
      GeoKeyDirectoryTag(1, 4) = NumberOfKeys;
      code = option.ProjLinearUnitSizeGeoKey;
      GeoKeyDirectoryTag = [GeoKeyDirectoryTag; [3077 34736 1 GeoDoubleOffset]];
      GeoDoubleParamsTag = [GeoDoubleParamsTag; code];
      GeoDoubleOffset = GeoDoubleOffset + 1;
   end

   % set ProjStdParallel1GeoKey, 3078
   if isfield(option, 'ProjStdParallel1GeoKey')
      NumberOfKeys = NumberOfKeys + 1;
      GeoKeyDirectoryTag(1, 4) = NumberOfKeys;
      code = option.ProjStdParallel1GeoKey;
      GeoKeyDirectoryTag = [GeoKeyDirectoryTag; [3078 34736 1 GeoDoubleOffset]];
      GeoDoubleParamsTag = [GeoDoubleParamsTag; code];
      GeoDoubleOffset = GeoDoubleOffset + 1;
   end

   % set ProjStdParallel2GeoKey, 3079
   if isfield(option, 'ProjStdParallel2GeoKey')
      NumberOfKeys = NumberOfKeys + 1;
      GeoKeyDirectoryTag(1, 4) = NumberOfKeys;
      code = option.ProjStdParallel2GeoKey;
      GeoKeyDirectoryTag = [GeoKeyDirectoryTag; [3079 34736 1 GeoDoubleOffset]];
      GeoDoubleParamsTag = [GeoDoubleParamsTag; code];
      GeoDoubleOffset = GeoDoubleOffset + 1;
   end

   % set ProjNatOriginLongGeoKey, 3080
   if isfield(option, 'ProjNatOriginLongGeoKey')
      NumberOfKeys = NumberOfKeys + 1;
      GeoKeyDirectoryTag(1, 4) = NumberOfKeys;
      code = option.ProjNatOriginLongGeoKey;
      GeoKeyDirectoryTag = [GeoKeyDirectoryTag; [3080 34736 1 GeoDoubleOffset]];
      GeoDoubleParamsTag = [GeoDoubleParamsTag; code];
      GeoDoubleOffset = GeoDoubleOffset + 1;
   end

   % set ProjNatOriginLatGeoKey, 3081
   if isfield(option, 'ProjNatOriginLatGeoKey')
      NumberOfKeys = NumberOfKeys + 1;
      GeoKeyDirectoryTag(1, 4) = NumberOfKeys;
      code = option.ProjNatOriginLatGeoKey;
      GeoKeyDirectoryTag = [GeoKeyDirectoryTag; [3081 34736 1 GeoDoubleOffset]];
      GeoDoubleParamsTag = [GeoDoubleParamsTag; code];
      GeoDoubleOffset = GeoDoubleOffset + 1;
   end

   % set ProjFalseEastingGeoKey, 3082
   if isfield(option, 'ProjFalseEastingGeoKey')
      NumberOfKeys = NumberOfKeys + 1;
      GeoKeyDirectoryTag(1, 4) = NumberOfKeys;
      code = option.ProjFalseEastingGeoKey;
      GeoKeyDirectoryTag = [GeoKeyDirectoryTag; [3082 34736 1 GeoDoubleOffset]];
      GeoDoubleParamsTag = [GeoDoubleParamsTag; code];
      GeoDoubleOffset = GeoDoubleOffset + 1;
   end

   % set ProjFalseNorthingGeoKey, 3083
   if isfield(option, 'ProjFalseNorthingGeoKey')
      NumberOfKeys = NumberOfKeys + 1;
      GeoKeyDirectoryTag(1, 4) = NumberOfKeys;
      code = option.ProjFalseNorthingGeoKey;
      GeoKeyDirectoryTag = [GeoKeyDirectoryTag; [3083 34736 1 GeoDoubleOffset]];
      GeoDoubleParamsTag = [GeoDoubleParamsTag; code];
      GeoDoubleOffset = GeoDoubleOffset + 1;
   end

   % set ProjFalseOriginLongGeoKey, 3084
   if isfield(option, 'ProjFalseOriginLongGeoKey')
      NumberOfKeys = NumberOfKeys + 1;
      GeoKeyDirectoryTag(1, 4) = NumberOfKeys;
      code = option.ProjFalseOriginLongGeoKey;
      GeoKeyDirectoryTag = [GeoKeyDirectoryTag; [3084 34736 1 GeoDoubleOffset]];
      GeoDoubleParamsTag = [GeoDoubleParamsTag; code];
      GeoDoubleOffset = GeoDoubleOffset + 1;
   end

   % set ProjFalseOriginLatGeoKey, 3085
   if isfield(option, 'ProjFalseOriginLatGeoKey')
      NumberOfKeys = NumberOfKeys + 1;
      GeoKeyDirectoryTag(1, 4) = NumberOfKeys;
      code = option.ProjFalseOriginLatGeoKey;
      GeoKeyDirectoryTag = [GeoKeyDirectoryTag; [3085 34736 1 GeoDoubleOffset]];
      GeoDoubleParamsTag = [GeoDoubleParamsTag; code];
      GeoDoubleOffset = GeoDoubleOffset + 1;
   end

   % set ProjFalseOriginEastingGeoKey, 3086
   if isfield(option, 'ProjFalseOriginEastingGeoKey')
      NumberOfKeys = NumberOfKeys + 1;
      GeoKeyDirectoryTag(1, 4) = NumberOfKeys;
      code = option.ProjFalseOriginEastingGeoKey;
      GeoKeyDirectoryTag = [GeoKeyDirectoryTag; [3086 34736 1 GeoDoubleOffset]];
      GeoDoubleParamsTag = [GeoDoubleParamsTag; code];
      GeoDoubleOffset = GeoDoubleOffset + 1;
   end

   % set ProjFalseOriginNorthingGeoKey, 3087
   if isfield(option, 'ProjFalseOriginNorthingGeoKey')
      NumberOfKeys = NumberOfKeys + 1;
      GeoKeyDirectoryTag(1, 4) = NumberOfKeys;
      code = option.ProjFalseOriginNorthingGeoKey;
      GeoKeyDirectoryTag = [GeoKeyDirectoryTag; [3087 34736 1 GeoDoubleOffset]];
      GeoDoubleParamsTag = [GeoDoubleParamsTag; code];
      GeoDoubleOffset = GeoDoubleOffset + 1;
   end

   % set ProjCenterLongGeoKey, 3088
   if isfield(option, 'ProjCenterLongGeoKey')
      NumberOfKeys = NumberOfKeys + 1;
      GeoKeyDirectoryTag(1, 4) = NumberOfKeys;
      code = option.ProjCenterLongGeoKey;
      GeoKeyDirectoryTag = [GeoKeyDirectoryTag; [3088 34736 1 GeoDoubleOffset]];
      GeoDoubleParamsTag = [GeoDoubleParamsTag; code];
      GeoDoubleOffset = GeoDoubleOffset + 1;
   end

   % set ProjCenterLatGeoKey, 3089
   if isfield(option, 'ProjCenterLatGeoKey')
      NumberOfKeys = NumberOfKeys + 1;
      GeoKeyDirectoryTag(1, 4) = NumberOfKeys;
      code = option.ProjCenterLatGeoKey;
      GeoKeyDirectoryTag = [GeoKeyDirectoryTag; [3089 34736 1 GeoDoubleOffset]];
      GeoDoubleParamsTag = [GeoDoubleParamsTag; code];
      GeoDoubleOffset = GeoDoubleOffset + 1;
   end

   % set ProjCenterEastingGeoKey, 3090
   if isfield(option, 'ProjCenterEastingGeoKey')
      NumberOfKeys = NumberOfKeys + 1;
      GeoKeyDirectoryTag(1, 4) = NumberOfKeys;
      code = option.ProjCenterEastingGeoKey;
      GeoKeyDirectoryTag = [GeoKeyDirectoryTag; [3090 34736 1 GeoDoubleOffset]];
      GeoDoubleParamsTag = [GeoDoubleParamsTag; code];
      GeoDoubleOffset = GeoDoubleOffset + 1;
   end

   % set ProjCenterOriginNorthingGeoKey, 3091
   if isfield(option, 'ProjCenterOriginNorthingGeoKey')
      NumberOfKeys = NumberOfKeys + 1;
      GeoKeyDirectoryTag(1, 4) = NumberOfKeys;
      code = option.ProjCenterOriginNorthingGeoKey;
      GeoKeyDirectoryTag = [GeoKeyDirectoryTag; [3091 34736 1 GeoDoubleOffset]];
      GeoDoubleParamsTag = [GeoDoubleParamsTag; code];
      GeoDoubleOffset = GeoDoubleOffset + 1;
   end

   % set ProjScaleAtNatOriginGeoKey, 3092
   if isfield(option, 'ProjScaleAtNatOriginGeoKey')
      NumberOfKeys = NumberOfKeys + 1;
      GeoKeyDirectoryTag(1, 4) = NumberOfKeys;
      code = option.ProjScaleAtNatOriginGeoKey;
      GeoKeyDirectoryTag = [GeoKeyDirectoryTag; [3092 34736 1 GeoDoubleOffset]];
      GeoDoubleParamsTag = [GeoDoubleParamsTag; code];
      GeoDoubleOffset = GeoDoubleOffset + 1;
   end

   % set ProjScaleAtCenterGeoKey, 3093
   if isfield(option, 'ProjScaleAtCenterGeoKey')
      NumberOfKeys = NumberOfKeys + 1;
      GeoKeyDirectoryTag(1, 4) = NumberOfKeys;
      code = option.ProjScaleAtCenterGeoKey;
      GeoKeyDirectoryTag = [GeoKeyDirectoryTag; [3093 34736 1 GeoDoubleOffset]];
      GeoDoubleParamsTag = [GeoDoubleParamsTag; code];
      GeoDoubleOffset = GeoDoubleOffset + 1;
   end

   % set ProjAzimuthAngleGeoKey, 3094
   if isfield(option, 'ProjAzimuthAngleGeoKey')
      NumberOfKeys = NumberOfKeys + 1;
      GeoKeyDirectoryTag(1, 4) = NumberOfKeys;
      code = option.ProjAzimuthAngleGeoKey;
      GeoKeyDirectoryTag = [GeoKeyDirectoryTag; [3094 34736 1 GeoDoubleOffset]];
      GeoDoubleParamsTag = [GeoDoubleParamsTag; code];
      GeoDoubleOffset = GeoDoubleOffset + 1;
   end

   % set ProjStraightVertPoleLongGeoKey, 3095
   if isfield(option, 'ProjStraightVertPoleLongGeoKey')
      NumberOfKeys = NumberOfKeys + 1;
      GeoKeyDirectoryTag(1, 4) = NumberOfKeys;
      code = option.ProjStraightVertPoleLongGeoKey;
      GeoKeyDirectoryTag = [GeoKeyDirectoryTag; [3095 34736 1 GeoDoubleOffset]];
      GeoDoubleParamsTag = [GeoDoubleParamsTag; code];
      GeoDoubleOffset = GeoDoubleOffset + 1;
   end

   % set VerticalCSTypeGeoKey, 4096
   if isfield(option, 'VerticalCSTypeGeoKey')
      NumberOfKeys = NumberOfKeys + 1;
      GeoKeyDirectoryTag(1, 4) = NumberOfKeys;
      code = option.VerticalCSTypeGeoKey;
      GeoKeyDirectoryTag = [GeoKeyDirectoryTag; [4096 0 1 code]];
   end

   % set VerticalCitationGeoKey, 4097
   if isfield(option, 'VerticalCitationGeoKey')
      NumberOfKeys = NumberOfKeys + 1;
      GeoKeyDirectoryTag(1, 4) = NumberOfKeys;
      code = [option.VerticalCitationGeoKey(:); '|']; cnt = length(code);
      GeoKeyDirectoryTag = [GeoKeyDirectoryTag; [4097 34737 cnt GeoAsciiOffset]];
      GeoAsciiParamsTag = [GeoAsciiParamsTag; code];
      GeoAsciiOffset = GeoAsciiOffset + cnt;
   end

   % set VerticalUnitsGeoKey, 4099
   if isfield(option, 'VerticalUnitsGeoKey')
      NumberOfKeys = NumberOfKeys + 1;
      GeoKeyDirectoryTag(1, 4) = NumberOfKeys;
      code = option.VerticalUnitsGeoKey;
      GeoKeyDirectoryTag = [GeoKeyDirectoryTag; [4099 0 1 code]];
   end

   % set ifd.ImageWidth, 256
   ifd.ImageWidth = size(image, 2);

   % set ifd.ImageLength, 257
   ifd.ImageLength = size(image, 1);

   % set ifd.BitsPerSample, 258
   ifd.BitsPerSample = abs(bit_depth)*ones(size(image,3),1);

   % set ifd.Compression, 259
   ifd.Compression = 1;				% Uncompressed

   % set ifd.PhotometricInterpretation, 262
   if isfield(option, 'ColorMap')
      ifd.PhotometricInterpretation = 3;	% Palette color
   else
      ifd.PhotometricInterpretation = 1;	% BlackIsZero
   end

   % set ifd.StripOffsets, 273
   ifd.StripOffsets = ones(ifd.ImageLength,1) * ifd.ImageWidth ...
	* abs(bit_depth)/8 * size(image,3);	% Num of bytes in a strip
   ifd.StripOffsets = cumsum(ifd.StripOffsets);
   ifd.StripOffsets(end) = [];
   ifd.StripOffsets = [0; ifd.StripOffsets];
   ifd.StripOffsets = ifd.StripOffsets + 8;	% 8 bytes before 1st strip

   % set ifd.SamplesPerPixel, 277
   ifd.SamplesPerPixel = size(image, 3);	% Gray level intensity

   % set ifd.RowsPerStrip, 278
   ifd.RowsPerStrip = 1;			% rows per strip

   % set ifd.StripByteCounts, 279
   ifd.StripByteCounts = ones(ifd.ImageLength,1) * ifd.ImageWidth ...
	* abs(bit_depth)/8 * size(image,3);	% Num of bytes in a strip

   % set ifd.XResolution, 282
   ifd.XResolution = 96;			% 96 dpi
   % set ifd.YResolution, 283
   ifd.YResolution = 96;			% 96 dpi

   % set ifd.PlanarConfiguration, 284
   ifd.PlanarConfiguration = 1;			% Chunky

   % set ifd.ResolutionUnit, 296
   ifd.ResolutionUnit = 2;			% Inch

   % set ifd.ColorMap, 320
   if isfield(option, 'ColorMap')
      ifd.ColorMap = option.ColorMap(:);
   end

   % set ifd.ExtraSamples, 338
   if ndims(image) > 2
      ifd.ExtraSamples = zeros((size(image,3)-1),1);
   end

   % set ifd.SampleFormat, 339
   switch bit_depth
   case {1, -1}
      ifd.SampleFormat = 4*ones(size(image,3),1);	% Undefined
   case {8, -16}
      ifd.SampleFormat = ones(size(image,3),1);		% Unsigned integer
   case {16, -32}
      ifd.SampleFormat = 2*ones(size(image,3),1);	% Signed integer
   case 32
      ifd.SampleFormat = 3*ones(size(image,3),1);	% Floating point
   end

   if isfield(ifd, 'ColorMap') && bit_depth == 16
      ifd.SampleFormat = ones(size(image,3),1);		% Unsigned integer
   end

   if ~isempty(bbox)

      % set gifd.ModelPixelScaleTag, 33550
      if gifd.GTRasterTypeGeoKey == 1
         Sx = (bbox(2) - bbox(1)) / ifd.ImageWidth;
         Sy = (bbox(4) - bbox(3)) / ifd.ImageLength;
      elseif gifd.GTRasterTypeGeoKey == 2
         Sx = (bbox(2) - bbox(1)) / (ifd.ImageWidth-1);
         Sy = (bbox(4) - bbox(3)) / (ifd.ImageLength-1);
      else
         error('Incorrect option.GTRasterTypeGeoKey value');
      end

      gifd.ModelPixelScaleTag = [Sx Sy 0];

      % set gifd.ModelTiepointTag, 33922
      RasterTiepoint = [0 0 0];
      ModelTiepoint = [bbox(1) bbox(4) 0];
      gifd.ModelTiepointTag = [RasterTiepoint ModelTiepoint];
   else
      if isfield(option, 'ModelPixelScaleTag')
         gifd.ModelPixelScaleTag = option.ModelPixelScaleTag(:);
      else
         error('ModelPixelScale is required when your Model Type is not Geographic.');
      end

      if isfield(option, 'ModelTiepointTag')
         gifd.ModelTiepointTag = option.ModelTiepointTag(:);
      end

      if isfield(option, 'ModelTransformationTag')
         gifd.ModelTransformationTag = option.ModelTransformationTag(:);
      end

      if ~isfield(option, 'ModelTiepointTag') && ~isfield(option, 'ModelTransformationTag')
         error('Either ModelTiepoint or ModelTransformation is required when your Model Type is not Geographic.');
      end
   end

   % write TIFF info
   num_entry = 16;				% 16 basic entries

   if isfield(ifd, 'ColorMap')
      num_entry = num_entry + 1;
   end

   if isfield(ifd, 'ExtraSamples')
      num_entry = num_entry + 1;
   end

   if isfield(gifd, 'ModelPixelScaleTag')
      num_entry = num_entry + 1;
   end

   if isfield(gifd, 'ModelTiepointTag')
      num_entry = num_entry + 1;
   end

   if isfield(gifd, 'ModelTransformationTag')
      num_entry = num_entry + 1;
   end

   if GeoAsciiOffset > 0
      num_entry = num_entry + 1;
   end

   if GeoDoubleOffset > 0
      num_entry = num_entry + 1;
   end

   ifd_end = 8 + sum(ifd.StripByteCounts) + 2 + num_entry*12 + 4;
   fid = fopen(filename, 'wb', 'ieee-le');

   fwrite(fid, 'II', 'uint8');			% Little-endian
   fwrite(fid, 42, 'uint16');			% TIFF signature
   fwrite(fid, (8 + sum(ifd.StripByteCounts)), 'uint32');	% IFD offset

   % write data strip
   image = squeeze(permute(image, [3 2 1]));

   switch abs(bit_depth)
   case 1
      image = bit2byte(image);
      fwrite(fid, image, 'uint8');
   case 8
      fwrite(fid, image, 'uint8');
   case 16
      if isfield(ifd, 'ColorMap') || bit_depth == -16
         fwrite(fid, image, 'uint16');
      else
         fwrite(fid, image, 'int16');
      end
   case 32
      if bit_depth == -32
         fwrite(fid, image, 'int32');
      else
         fwrite(fid, image, 'single');
      end
   end

   fwrite(fid, num_entry, 'uint16');

   %  Entry 1
   fwrite(fid, 256, 'uint16');
   fwrite(fid, 3, 'uint16');		% uint16
   fwrite(fid, 1, 'uint32');		% count
   fwrite(fid, ifd.ImageWidth, 'uint16');
   fwrite(fid, 0, 'uint16');

   %  Entry 2
   fwrite(fid, 257, 'uint16');
   fwrite(fid, 3, 'uint16');		% uint16
   fwrite(fid, 1, 'uint32');		% count
   fwrite(fid, ifd.ImageLength, 'uint16');
   fwrite(fid, 0, 'uint16');

   %  Entry 3
   fwrite(fid, 258, 'uint16');
   fwrite(fid, 3, 'uint16');		% uint16
   fwrite(fid, length(ifd.BitsPerSample), 'uint32');	% count

   if length(ifd.BitsPerSample) > 2
      fwrite(fid, ifd_end, 'uint32');
      ifd_end = ifd_end + length(ifd.BitsPerSample)*2;
   else
      fwrite(fid, ifd.BitsPerSample(1), 'uint16');

      if length(ifd.BitsPerSample) > 1
         fwrite(fid, ifd.BitsPerSample(2), 'uint16');
      else
         fwrite(fid, 0, 'uint16');
      end
   end

   %  Entry 4
    fwrite(fid, 259, 'uint16');
   fwrite(fid, 3, 'uint16');		% uint16
   fwrite(fid, 1, 'uint32');		% count
   fwrite(fid, ifd.Compression, 'uint16');
   fwrite(fid, 0, 'uint16');

   %  Entry 5
   fwrite(fid, 262, 'uint16');
   fwrite(fid, 3, 'uint16');		% uint16
   fwrite(fid, 1, 'uint32');		% count
   fwrite(fid, ifd.PhotometricInterpretation, 'uint16');
   fwrite(fid, 0, 'uint16');

   %  Entry 6
   fwrite(fid, 273, 'uint16');		% StripOffsets
   fwrite(fid, 4, 'uint16');		% uint32
   fwrite(fid, ifd.ImageLength, 'uint32');	% count
   fwrite(fid, ifd_end, 'uint32');
   ifd_end = ifd_end + ifd.ImageLength*4;

   %  Entry 7
   fwrite(fid, 274, 'uint16');
   fwrite(fid, 3, 'uint16');		% uint16
   fwrite(fid, 1, 'uint32');		% count
   fwrite(fid, ifd.Orientation, 'uint16');
   fwrite(fid, 0, 'uint16');

   %  Entry 8
   fwrite(fid, 277, 'uint16');
   fwrite(fid, 3, 'uint16');		% uint16
   fwrite(fid, 1, 'uint32');		% count
   fwrite(fid, ifd.SamplesPerPixel, 'uint16');
   fwrite(fid, 0, 'uint16');

   %  Entry 9
   fwrite(fid, 278, 'uint16');
   fwrite(fid, 3, 'uint16');		% uint16
   fwrite(fid, 1, 'uint32');		% count
   fwrite(fid, ifd.RowsPerStrip, 'uint16');
   fwrite(fid, 0, 'uint16');

   %  Entry 10
   fwrite(fid, 279, 'uint16');		% StripByteCounts
   fwrite(fid, 4, 'uint16');		% uint32
   fwrite(fid, ifd.ImageLength, 'uint32');	% count
   fwrite(fid, ifd_end, 'uint32');
   ifd_end = ifd_end + ifd.ImageLength*4;

   %  Entry 11
   fwrite(fid, 282, 'uint16');		% XResolution
   fwrite(fid, 5, 'uint16');		% uint16
   fwrite(fid, 1, 'uint32');		% count
   fwrite(fid, ifd_end, 'uint32');
   ifd_end = ifd_end + 2*4;

   %  Entry 12
   fwrite(fid, 283, 'uint16');		% YResolution
   fwrite(fid, 5, 'uint16');		% uint16
   fwrite(fid, 1, 'uint32');		% count
   fwrite(fid, ifd_end, 'uint32');
   ifd_end = ifd_end + 2*4;

   %  Entry 13
   fwrite(fid, 284, 'uint16');
   fwrite(fid, 3, 'uint16');		% uint16
   fwrite(fid, 1, 'uint32');		% count
   fwrite(fid, ifd.PlanarConfiguration, 'uint16');
   fwrite(fid, 0, 'uint16');

   %  Entry 14
   fwrite(fid, 296, 'uint16');
   fwrite(fid, 3, 'uint16');		% uint16
   fwrite(fid, 1, 'uint32');		% count
   fwrite(fid, ifd.ResolutionUnit, 'uint16');
   fwrite(fid, 0, 'uint16');

   if isfield(option, 'ColorMap')
      fwrite(fid, 320, 'uint16');
      fwrite(fid, 3, 'uint16');		% uint16
      fwrite(fid, length(ifd.ColorMap), 'uint32');	% count
      fwrite(fid, ifd_end, 'uint32');
      ifd_end = ifd_end + length(ifd.ColorMap)*2;
   end

   if isfield(ifd, 'ExtraSamples')
      fwrite(fid, 338, 'uint16');
      fwrite(fid, 3, 'uint16');		% uint16
      fwrite(fid, length(ifd.ExtraSamples), 'uint32');	% count

      if length(ifd.ExtraSamples) > 2
         fwrite(fid, ifd_end, 'uint32');
         ifd_end = ifd_end + length(ifd.ExtraSamples)*2;
      else
         fwrite(fid, ifd.ExtraSamples(1), 'uint16');

         if length(ifd.ExtraSamples) > 1
            fwrite(fid, ifd.ExtraSamples(2), 'uint16');
         else
            fwrite(fid, 0, 'uint16');
         end
      end
   end

   %  Entry 15
   fwrite(fid, 339, 'uint16');
   fwrite(fid, 3, 'uint16');		% uint16
   fwrite(fid, length(ifd.SampleFormat), 'uint32');	% count

   if length(ifd.SampleFormat) > 2
      fwrite(fid, ifd_end, 'uint32');
      ifd_end = ifd_end + length(ifd.SampleFormat)*2;
   else
      fwrite(fid, ifd.SampleFormat(1), 'uint16');

      if length(ifd.SampleFormat) > 1
         fwrite(fid, ifd.SampleFormat(2), 'uint16');
      else
         fwrite(fid, 0, 'uint16');
      end
   end

   if isfield(gifd, 'ModelPixelScaleTag')
      fwrite(fid, 33550, 'uint16');	% ModelPixelScaleTag
      fwrite(fid, 12, 'uint16');	% double
      fwrite(fid, 3, 'uint32');		% count
      fwrite(fid, ifd_end, 'uint32');
      ifd_end = ifd_end + 3*8;
   end

   if isfield(gifd, 'ModelTiepointTag')
      fwrite(fid, 33922, 'uint16');	% ModelTiepointTag
      fwrite(fid, 12, 'uint16');	% double
      fwrite(fid, length(gifd.ModelTiepointTag), 'uint32');	% count
      fwrite(fid, ifd_end, 'uint32');
      ifd_end = ifd_end + length(gifd.ModelTiepointTag)*8;
   end

   if isfield(gifd, 'ModelTransformationTag')
      fwrite(fid, 34264, 'uint16');	% ModelTransformationTag
      fwrite(fid, 12, 'uint16');	% double
      fwrite(fid, 16, 'uint32');	% count
      fwrite(fid, ifd_end, 'uint32');
      ifd_end = ifd_end + 16*8;
   end

   %  Entry 16
   fwrite(fid, 34735, 'uint16');	% GeoKeyDirectoryTag
   fwrite(fid, 3, 'uint16');		% uint16
   fwrite(fid, numel(GeoKeyDirectoryTag), 'uint32');	% count
   fwrite(fid, ifd_end, 'uint32');
   ifd_end = ifd_end + numel(GeoKeyDirectoryTag)*2;

   if GeoDoubleOffset > 0
      fwrite(fid, 34736, 'uint16');	% GeoDoubleParamsTag
      fwrite(fid, 12, 'uint16');	% double
      fwrite(fid, length(GeoDoubleParamsTag), 'uint32');	% count
      fwrite(fid, ifd_end, 'uint32');
      ifd_end = ifd_end + length(GeoDoubleParamsTag)*8;
   end

   if GeoAsciiOffset > 0
      GeoAsciiParamsTag = [GeoAsciiParamsTag; 0];	% TIFF6 requres NUL ending
      fwrite(fid, 34737, 'uint16');	% GeoAsciiParamsTag
      fwrite(fid, 2, 'uint16');		% ascii
      fwrite(fid, length(GeoAsciiParamsTag), 'uint32');		% count
      fwrite(fid, ifd_end, 'uint32');
   end

   %  IFD is terminated with 4-byte offset to the next IFD, or 0 if there are none.
   fwrite(fid, 0, 'uint32');

   if length(ifd.BitsPerSample) > 2
      fwrite(fid, ifd.BitsPerSample, 'uint16');		% 258
   end

   fwrite(fid, ifd.StripOffsets, 'uint32');		% 273
   fwrite(fid, ifd.StripByteCounts, 'uint32');		% 279
   fwrite(fid, ifd.XResolution, 'uint32');		% 282
   fwrite(fid, 1, 'uint32');				% 282
   fwrite(fid, ifd.YResolution, 'uint32');		% 283
   fwrite(fid, 1, 'uint32');				% 283

   if isfield(option, 'ColorMap')
      fwrite(fid, ifd.ColorMap, 'uint16');		% 320
   end

   if isfield(ifd, 'ExtraSamples') && length(ifd.ExtraSamples) > 2
      fwrite(fid, ifd.ExtraSamples, 'uint16');		% 338
   end

   if length(ifd.SampleFormat) > 2
      fwrite(fid, ifd.SampleFormat, 'uint16');		% 339
   end

   fwrite(fid, gifd.ModelPixelScaleTag, 'double');	% 33550
   fwrite(fid, gifd.ModelTiepointTag, 'double');	% 33922
   GeoKeyDirectoryTag = GeoKeyDirectoryTag';
   fwrite(fid, GeoKeyDirectoryTag, 'uint16');		% 34735

   if GeoDoubleOffset > 0
      fwrite(fid, GeoDoubleParamsTag, 'double');	% 34736
   end

   if GeoAsciiOffset > 0
      fwrite(fid, GeoAsciiParamsTag, 'uint8');		% 34737
   end

   fclose(fid);

   return;					

function byte = bit2byte(bit)

   byte = reshape(bit, [8, length(bit(:))/8])';
   byte(byte>0) = 1;
   byte(byte<1) = 0;
   byte = double(byte)*[128 64 32 16 8 4 2 1]';

return;	