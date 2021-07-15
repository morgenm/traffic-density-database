%{
Extract Cells:
This program extracts the risk (event rate) for given latitude and longitude
ranges. It outputs the risks for each cell (in MIT's cell format) to a
CSV file. It additionally outputs the latitude and longitude bounds for
each cell to a separate CSV file, so it can be converted to our cell
format. This program is to be used with the MIT TrafficDensityDatabase.
%}

%{ 
This program extracts risk data from the low altitude MIT TrafficDensity
database. The code for this is pulled from MIT's example code here:
https://github.com/mit-ll/traffic-density-database/blob/main/Examples/simpleExampleLowAltitude.m
%}
td = TrafficDensity('joblowh.dat');

% Specify low-altitude tables
td.filenames.cell = 'celllowh.mat';
td.filenames.cellCoverage = 'cellcoveragelowh.mat';
td.filenames.cellAirspace = 'cellAirspacelowh.mat';

td = td.LoadData;

% The area we need to find risk for
td.area.LongitudeLimit = [-90.300, -89.900];
td.area.LatitudeLimit = [32.100, 32.600];

td.plotresults = false; % Save plotting for below

td = td.run;

td.processNoncoop = true; % Add plotting of noncooperatives (it is only the plotting function that considers noncooperatives)

% Plot density with all options (see plot method documentation for other parameters)
% The outdata and summarize outputs provide user access to the estimated parameters
[outdata,summarize] = td.plot('plottype','rate',...
    'plotscale','linear',...
    'plotvalue','avg',...
    'plotcombined',true,...
    'plottitleadd',' (exampleAreaPlot)',...
    'returnonlydata',false,...
    'summarizearea',true);

% Write the aggregate risk for 100-500m to CSV file
% outdata(4) is aggregate and the the first rateavg is for 100-500m
riskAggregateLowAlt = outdata(4).rateavg(:,:,1);
riskAggregateLowAlt = flip(riskAggregateLowAlt, 1);
writematrix(riskAggregateLowAlt, "riskAggregateLowAlt.csv");

% Find all bounding coordinates for each cell in the risk area
numCellsLong = td.cellLonLim(2)-td.cellLonLim(1);
numCellsLat = td.cellLatLim(2)-td.cellLatLim(1);
cellBound = zeros(numCellsLat, numCellsLong);
cellBound = string(cellBound); % Create string matrix full of zeros
for long = td.cellLonLim(1):td.cellLonLim(2)
    for lat = td.cellLatLim(1):td.cellLatLim(2)
        % Get the left, right, bottom and top bounds for the current cell
        lLong = td.cellLonCutpoints(long);
        rLong = td.cellLonCutpoints(long+1);
        bLat = td.cellLatCutpoints(lat);
        tLat = td.cellLatCutpoints(lat+1);
        
        % Create string of bounding coordinates to output to CSV
        blBound = [lLong, bLat];
        trBound = [rLong, tLat];
        bounds = [blBound, trBound];
        boundStrArr = string(bounds);
        boundStr = "(" + boundStrArr(1) + ", " + boundStrArr(2) + ")";
        boundStr = boundStr + ", (" + boundStrArr(3) + ", " + boundStrArr(4) + ")";
        
        % Add it to the matrix
        currLongCell = (long - td.cellLonLim(1))+1;
        currLatCell = (lat - td.cellLatLim(1))+1;
        cellBound(currLatCell, currLongCell) = boundStr;
    end
end

% Write cell bound matrix to CSV
cellBound = flip(cellBound, 1);
writematrix(cellBound, "cellBound.csv");

% Get MIT cell numbers
mitCells = zeros(numCellsLat, numCellsLong);
mitCells = string(mitCells);
for long = td.cellLonLim(1):td.cellLonLim(2)
    for lat = td.cellLatLim(1):td.cellLatLim(2)
        currLongIndex = (long - td.cellLonLim(1))+1;
        currLatIndex = (lat - td.cellLatLim(1))+1;
        longStr = string(long);
        latStr = string(lat);
        longLatStr = "(" + longStr + "," + latStr + ")";
        mitCells(currLatIndex, currLongIndex) = longLatStr;
    end
end

mitCells = flip(mitCells, 1);
writematrix(mitCells, "mitCellIndices.csv")
