function im = import_tif (filename)
% IMPORT_TIF Imports a .tif file into the workspace as an i-by-j-by-f
%   double array
%
%   im = import_tif uses GUI to select the file
%   im = import_tif(filename) loads .tif file from filename
%
%   Examples
%   --------
%   >> Fibronectin = import_tif;
%
%   See also IMPORT_TIF_SEQUENCE

if nargin < 1
    [filename, pathname] = uigetfile({'*.*'});
    filename = strcat(pathname,filename);
    disp(filename)
end

info = imfinfo(filename);
frames = size(info,1);
frames = max(frames, length(info.BitsPerSample));
sizei = info.Height;
sizej = info.Width;
im = nan(sizei,sizej,frames);

% if nargin < 2
%     [filename, pathname] = uigetfile({'*.*'});
% end


% for f=1:frames
%     im(:,:,f) = imread(filename,f);
% end

% Optional: plot
% image(x,y,im, 'CDataMapping', 'scaled');
% colormap(gray);
% caxis([min(im(:)) prctile(im(:),99)]);
% axis equal;

im = imread(filename);
im = uint16(im);

end
