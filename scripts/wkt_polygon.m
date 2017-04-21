function [ wkt ] = wkt_polygon( xxyy )
% Generates a num2string for a WKT polygon, to open as a new delimited text
% layer in QGIS
% 
% Input 
%   xxyy: Vector containing the min and max X and Y easting/northing
%         coordinates of the form [xMin xMax yMin yMax]
%        
%  Output
%   wkt: string defining square polygon extent for CSV 



wkt = ['POLYGON((', num2str(xxyy(1)), ' ', num2str(xxyy(3)),', ',...
                    num2str(xxyy(2)), ' ', num2str(xxyy(3)),', ',...
                    num2str(xxyy(2)), ' ', num2str(xxyy(4)),', ',...
                    num2str(xxyy(1)), ' ', num2str(xxyy(4)),', ',...
                    num2str(xxyy(1)), ' ', num2str(xxyy(3)),'))'];


end

