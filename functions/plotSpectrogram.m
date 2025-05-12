function [figHandle] = plotSpectrogram(...
  data,timeAxis,velVec,dbRange,flipData,extra,figHandle)
% [figHandle] = plotSpectrogram(data,timeAxis,velVec,dbRange,flipData,extra)
%
% Plotting the spectrogram in dB-scale. This code can also be used for data from 
% multiple images, where each image is a matrix of data, i.e., image nbr k
% is stored as data(:,:,k).
%
% IN:
%   data     : matrix with the data to plot. In the case of multiple images,
%              store each image as a matrix, with the image number on the
%              3rd dimension (image k: data(:,:,k));
%   timeAxis : vector (matrix in the case of multiple images, where the 
%              timings for the images are stored as columns) for the 
%              timings of the spectrogram, i.e. the x-axis of the plot;
%   velVec   : y-axis of the plot;
%   dbRange  : the dB range of the plot, i.e., the scaling of the 
%              power density;
%   flipData : 1 - reverse the y-axis in the spectrogram
%              0 - no flip
%              Needed since some functions outputs the spectrum in the 
%              wrong direction;
%   extra    : These can be omitted:
%              .nbrOfLevels - nbr of color levels in plot
%              .figTitle - title of plot
%
% OUT: 
%  figHandle : figure handle


try 
  nbrOfLevels = extra.nbrOfLevels;
catch
  nbrOfLevels = 128;
end

if nargin < 7
  figHandle = figure();
end
%   nbrOfLevels = 512;

nbrOfFreqs = size(data,1);
noImages = size(data,3);

if noImages == 1
  %Make sure the time axis is a row vector
  timeAxis = timeAxis(:);
end

%Turn data into dB-scale
data_db = 10*log10( data );
%Find the maximum intensity over all images, to be used for scaling the intensity correctly
maxData_db = max(max(max( data_db)));

hold on
for ii = 1:noImages
  
  data_dbTmp = double(( (( data_db(:,:,ii) - maxData_db ) + dbRange )/dbRange*nbrOfLevels));
  
  if flipData
    image(timeAxis(:,ii),velVec,data_dbTmp(nbrOfFreqs:-1:1,:));
  else
    image(timeAxis(:,ii),velVec,data_dbTmp);
  end
  
end

hold off

try
  title(extra.figTitle, 'FontSize', 16);
catch
  %Do nothing
end
  
xlabel('Time [s]');
ylabel('Velocity [m/s]', 'FontSize', 14);
set(gca,'YDir','normal', 'FontSize', 14)

% Set tick label font size (numbers on the axes)
ax = gca;                                % Get current axes handle
set(ax, 'FontSize', 13);                 % Set font size of tick labels
axis([0 0.18 min(velVec) max(velVec)])
% axis([0 0.3 -0.22 0.22])

colormap(gray(nbrOfLevels))
if exist('boldify','file')
  boldify
end

if nargin < 7
  figPos = get(figHandle,'Position');
  set(figHandle,'Position', [figPos(1) figPos(2) figPos(3) figPos(4)/2])
end