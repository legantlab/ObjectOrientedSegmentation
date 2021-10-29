function outstruct = calculate_order_parameter (varargin)
% CALCULATE_ORDER_PARAMETER analyzes a fluorescent image of F-actin and
%   calculates the nematic order parameter from an image of actin, as well
%   as the area, coherence, aspect ratio, and local orientation
%
%   Outstruct = calculate_order_parameter(I);
%       Analyzes the image I for orientations, order parameter, and
%       coherence
%   Outstruct = calculate_order_parameter(I, N1, V1, ...) is the
%       same with some additional N/V name/value pairs
%
%   Optional name/value pairs
%   -------------------------
%   'nmpp' -- nanometers per pixel of the image (default is 110), this is
%       not needed unless the area or focal adhesions are needed
%   'threshold' -- for thresholding, default=0.04, increase if there are
%       too many pixels and decrease if there are too few
%   'coherence_cutoff' -- removes some pixels that have unreliable
%       orientation (default=0.08)
%   'grdsigma' -- default=1, sigma of a Gaussian for calculating
%       gradients of structure tensor (analogous to window size)
%   'blksigma' -- default=3, sigma of Gaussian for smoothing data
%   'weighted' -- 'yes' for weighted, 'no' for unweighted, default = 'yes'
%   'weights_nptophat' -- for weighted calculations a top-hat filter will
%       be run, this is the number of pixels in the top-hat (default=8)
%   'maxweightprctile' -- (default=0.92) Percentile of the maximum pixel
%       value for weighting, so a few bright pixels do not dominate
%   'weights_grd' -- default = 0.6, Gaussian blur gradient during weight
%   'plot' -- show some intermediate plots, good for debugging, options
%       are 'yes' and 'no', default is 'no', if 'yes' is chosen then it
%       will save an image 'calculate_order_parameter.bmp' in the current
%       directory
%   'crop' -- allows you to freehand draw the cell before the function,
%       options are 'yes' and 'no', default = 'no'
%   'cropmask' -- crops to a prexisted mask
%   'FA_analysis' -- 'yes' to do optional analysis on focal adhesions,
%       there are some other parameters than be changed associated with
%       this. This is experimental, use at your own risk, change
%       'plot' to 'yes' if you try it
%   
%   Output structure
%   ----------------
%   S -- the order parameter
%   S_radial -- order parameter in radial coordinates
%   area -- the cell area
%   aspect_ratio -- the aspect ratio (calculated from the area)
%   avg_coherence -- the average coherence of the structure tensor
%   orientation -- the average orientation of actin (xy coordinates)
%   com -- center of mass (ij coordinates)
%   Orientim -- image of all local orientations from structure tensor
%   Coherence -- image of coherence from structure tensor
%   Weights -- image of weights used in calculations
%   Mask -- image of a mask
%   parameters -- input parameters to the program, for bookkeeping
%   FA_results -- in case FA_analysis was on
%   
%   Examples
%   --------
%   >> Outstruct = calculate_order_parameter(image, 'plot', 'yes');
%   >> Outstruct = calculate_order_parameter(image, 'threshold', 0.02);
%
%   Notes
%   -----
%   Bryant Doss, Benoit Ladoux
%   Mechanobiology Institute, National University of Singapore
%   requires RIDGESEGMENT, RIDGEORIENT, NORMALISE, DERIVATIVE5, and 
%   GAUSSFILT by Peter Kovesi, and EXCISE
%   adapted from code by Mukund Gupta
%   Version date: 2019-08-23
%   Changelog: 2019-08-23 added gradient input during weighting step

p = inputParser;
p.addRequired('Input'); % image of actin
p.addOptional('nmpp', 110, @isnumeric); % nanometers per pixel
p.addOptional('plot', 'no', @ischar); % plot the results
p.addOptional('coherence_cutoff', 0.08, @isnumeric);
p.addOptional('threshold', 0.04, @isnumeric); % lower for more pixels in mask, raise for less
p.addOptional('grdsigma', 1); % for ridgeorient
p.addOptional('blksigma', 3); % for ridgeorient
p.addOptional('crop', 'no', @ischar); % for cropping a cell
p.addOptional('cropmask', [], @islogical); % if we already have a mask
p.addOptional('weighted', 'yes', @ischar);
p.addOptional('weights_nptophat', 8, @isnumeric);
p.addOptional('maxweightprctile', 0.92, @isnumeric);
p.addOptional('weights_grd', 0.6, @isnumeric);
% Options for focal adhesion analysis
p.addOptional('FA_analysis', 'no', @ischar);
p.addOptional('FA_threshold', 2, @isnumeric);
p.addOptional('FA_minsize', 0.25, @isnumeric); % units: um^2
p.addOptional('FA_mergesize', 0.50, @isnumeric); % units: um^2
p.parse(varargin{:});
Input = p.Results.Input;
nmpp = p.Results.nmpp;
doplot = p.Results.plot;
coherence_cutoff = p.Results.coherence_cutoff;
grdsigma = p.Results.grdsigma;
blksigma = p.Results.blksigma;
thresh = p.Results.threshold;
docrop = p.Results.crop;
cropmask = p.Results.cropmask;
doweighted = p.Results.weighted;
weights_nptophat = p.Results.weights_nptophat;
maxweightprctile = p.Results.maxweightprctile;
weights_grd = p.Results.weights_grd;
FA_analysis = p.Results.FA_analysis;
FA_threshold = p.Results.FA_threshold;
FA_minsize = p.Results.FA_minsize;
FA_mergesize = p.Results.FA_mergesize;
outstruct.parameters = p.Results;

% Convert image to double
im = double(Input);
% Save the original input image
input_im = im;

% OPTIONAL: crop the image
if strcmp(docrop, 'yes')
    figure(9)
    image(im,'Cdatamapping','scaled')
    axis image
    h = imfreehand;
    wait(h);
    cropmask(:,:) =createMask(h);
    im = im .* cropmask;
    cropdist(:,:) = bwdist(abs(cropmask-1));
    close(9)
elseif ~isempty(cropmask)
    im = im .* cropmask;
    cropdist(:,:) = bwdist(abs(cropmask-1));
end


% Segment and normalize the image (for masking)
disp('Ridge Segment')
blksze = fix(2*blksigma); % just use this, seems to work fine
[normim, mask] = ridgesegment(im, blksze, thresh);


% Calculate the local orientation of actin using local structure tensor
disp('Ridge Orient')
[orientim, ~, coherence] = ridgeorient(normim, grdsigma, blksigma, 0);


% Kill things on the edge if it was cropped
if strcmp(docrop, 'yes')
    for i=1:size(cropdist,1)
        for j=1:size(cropdist,2)
            if cropdist(i,j) < blksze*3
                mask(i,j) = 0;
                orientim(i,j) = nan;
                coherence(i,j) = nan;
            end
        end
    end
end


% Kill things on the edge of the frame, not trustworthy data
disp(['Removing ' num2str(blksze*3) ' pixels from edge'])
for i=1:size(im,1)
    for j=1:size(im,2)
        if i < blksze*3 || i > size(im,1)-blksze*3 || j < blksze*3 || j > size(im,2)-blksze*3
            mask(i,j) = 0;
            orientim(i,j) = nan;
            coherence(i,j) = nan;
        end
    end
end


% "Coherence" from the orientation image, and change 0s to nan
for i=1:size(orientim,1)
    for j=1:size(orientim,2)
        if orientim(i,j)==0 || mask(i,j)==0 || coherence(i,j) < coherence_cutoff
            orientim(i,j)=nan;
            coherence(i,j) = nan;
            mask(i,j) = 0;
        end
    end
end


% Select only largest connected pixel area and remove the others.
% Basically keeps only one cell in the center of the image
regions = regionprops(mask, 'PixelList');
len = 0;
for a=1:length(regions)
    len = max(size(regions(a).PixelList,1),len);
end
for a=1:length(regions)
    if size(regions(a).PixelList,1) ~= len
        for b=1:size(regions(a).PixelList,1)
            mask(regions(a).PixelList(b,2),regions(a).PixelList(b,1)) = 0;
            orientim(regions(a).PixelList(b,2),regions(a).PixelList(b,1)) = nan;
            coherence(regions(a).PixelList(b,2),regions(a).PixelList(b,1)) = nan;
        end
    end
end
outstruct.Orientim = orientim;
outstruct.Coherence = coherence;

% Calculate the area of the cell, in um^2
outstruct.area = sum(sum(imfill(mask, 'holes'))) .* (nmpp/1000).^2;

% Apect ratio calculation
AR_results=regionprops(imfill(mask,'holes'), 'MajorAxisLength', ...
        'MinorAxisLength','Centroid','Orientation');
if ~isempty(AR_results)
    outstruct.aspect_ratio = AR_results.MajorAxisLength./AR_results.MinorAxisLength;
else
    outstruct.aspect_ratio = nan;
end
outstruct.Mask = imfill(mask,'holes');

% Filter the input image, at this point "im" is only used for weighting the
% order parameter and for displaying the final image
% Make the weights for S, it is "values"
masknan = double(mask);
for i=1:size(mask,1)
    for j=1:size(mask,2)
        if mask(i,j) == 0; masknan(i,j) = nan; end
    end
end
if strcmp(doweighted, 'yes')
    disp('Performing top-hat filter for weighting of S...');
    % First do a quick Gaussian filter
    im = gaussfilt(im, weights_grd);
    % Then do a top-hat filter
    im = imtophat(im, strel('disk', weights_nptophat));
    values = excise(masknan(:).*im(:));
    values = values - min(values);
    % Set a maximum weight so a few bright pixels won't overpower the image
    maxval = prctile(values,maxweightprctile*100);
    for v=1:length(values); if values(v)>maxval; values(v)=maxval; end; end
    S_calc_weights = masknan .* im;
    S_calc_weights(S_calc_weights>maxval) = maxval;
    % Use this for unweighted
elseif strcmp(doweighted, 'no');
    values = excise(masknan(:));
    S_calc_weights = masknan;
end
outstruct.Weights = S_calc_weights;
outstruct.avg_coherence = nansum(coherence(:).*S_calc_weights(:))./nansum(S_calc_weights(:));
disp(['Average Coherence: ' num2str(outstruct.avg_coherence)]);


% Calculate the global order parameter
angles = excise(orientim(:));
St = zeros(1, 180); avg_angle = zeros(1, 180);
% Try different rotations (so that the average angle is about pi/2, this
% is actually maximizing S but the effect is effectively the same)
for is=1:1:180
    angles_new = angles + is*pi/180;
    for is2=1:length(angles_new)
        if angles_new(is2) > pi;
            angles_new(is2) = angles_new(is2) - pi;
        end
    end
    avg = nansum(angles_new.*values)./nansum(values); % weighted mean
    avg_angle(is) = avg;
    diff = abs(avg - angles_new);
    St(is) = nansum(cos(2.*diff).*values)./nansum(values); % weighted mean
end
[S_xy,i] = max(St);
outstruct.S = S_xy;
outstruct.orientation = avg_angle(i)-i*pi/180; % average angle
disp(['Nematic Order Parameter (xy): ' num2str(S_xy)]);


% Calculate radial order parameter
% First calculate the (weighted or unweighted) center of mass
cellmask = S_calc_weights;
[jmat, imat] = meshgrid(1:size(cellmask, 2), 1:size(cellmask, 1));
weightedI = imat .* cellmask;
weightedJ = jmat .* cellmask;
iCOM = nansum(weightedI(:)) / nansum(excise(cellmask(:)));
jCOM = nansum(weightedJ(:)) / nansum(excise(cellmask(:)));
selectCOM = 'no'; % change this hardcode to select center of mass
if strcmp(selectCOM, 'yes')
    figure(2010); 
    imshow(im,[]);
    [jCOM,iCOM]=ginput(1);
    close(2010);
end
% Calculate distance and angle of each pixel from the center of mass
distance_com = nan(size(orientim,1), size(orientim,2));
angle_com = nan(size(orientim,1), size(orientim,2));
radial_orient = nan(size(orientim,1), size(orientim,2));
for i=1:size(orientim,1)
    for j=1:size(orientim,2)
        if isnan(orientim(i,j))
            distance_com(i,j) = nan;
            angle_com(i,j) = nan;
            continue
        end
        distance_com(i,j) = norm([i-iCOM j-jCOM]);
        angle_com(i,j) = -atan2(j-jCOM, i-iCOM) + pi/2;
        radial_orient(i,j) = cos(2.*(orientim(i,j) - angle_com(i,j)));        
    end
end
S_radial = excise(radial_orient(:));
S_radial = nansum(S_radial .* values) ./ nansum(values);
outstruct.S_radial = S_radial;
% disp(['Nematic Order Parameter (radial): ' num2str(S_radial)]);
COM = [iCOM jCOM];
outstruct.com = COM;


%%% Statistics for focal adhesions, watershed algorithm
% Based on algorithm from: Zamir et al., J. Cell Sci. 1999, 1655-1669
if strcmp(FA_analysis, 'yes')
    disp('Performing optional watershed analysis on focal adhesions...');
    % Process input image
    im = input_im;
    im = imfilter(im, fspecial('gaussian', 5, .6));
    im = imtophat(im, strel('disk', 16));
    pix = FA_minsize / (nmpp * nmpp) * 1e6;  % .25um^2 cutoff
    size_tol = ceil(pix); % tolerance for focal adhesion size
    pix = FA_mergesize / (nmpp * nmpp) * 1e6;  % .50um^2 cutoff
    size_tol2 = ceil(pix);% tolerance for combining adhesions
    % Define an intensity threshold for performing the watershed
    threshold = min(min(blkproc(im, [ceil(size(im,1)/8) ceil(size(im,2)/8)], @(x) max(x(:)))));
    threshold = threshold * FA_threshold; % Intensity threshold for the watershed
    if strcmp(docrop,'yes'); im = im .* cropmask; end % ignore cropped pixels if needed
    disp(['Pixel thresholds: ' num2str(size_tol) ', ' num2str(size_tol2) ...
        ', intensity threshold: ' num2str(threshold)]);
    % We are already working with im which has been filtered
    % Sort the pixels so the brightest are first (descending)
    [values, indices] = sort(im(:), 'descend');
    shed = nan(size(im));
    shedmax = 1;
    shedsizes = zeros(1,nansum(nansum(values>threshold)));
    % Loop over all pixels in the image until we are below the threshold
    for iv=1:length(values)
        if values(iv) < threshold
            break
        end
        % Build up a neighborhood of nearby pixels to see which watershed
        % they are in
        ip = mod(indices(iv)-1, size(im,1)) + 1;
        jp = ceil(indices(iv) / size(im,1));
        % Neighbors are diagonal and top/bottom and left/right
        neighbors = [ip+1 jp; ip-1 jp; ip jp+1; ip jp-1; ip+1 jp+1; ...
                     ip-1 jp+1; ip+1 jp-1; ip-1 jp-1];
        shedneighbors = [];
        % Loop over the neighbors, see which shed they belong to
        for in=1:size(neighbors,1)
            if neighbors(in,1) < size(shed,1) && neighbors(in,1) > 0 && ...
                    neighbors(in,2) < size(shed,2) && neighbors(in,2) > 0
                shedneighbors = [shedneighbors shed(neighbors(in,1), neighbors(in,2))];
            end
        end
        shedneighbors = excise(shedneighbors); % remove NaN
        shedneighbors = unique(shedneighbors); % remove duplicates
        if isempty(shedneighbors) % if no neighbors, create a new shed
            shed(ip,jp) = shedmax;
            shedsizes(shedmax) = 1;
            shedmax = shedmax + 1;
        elseif length(shedneighbors) == 1 % if only one neighbor, append
            shed(ip,jp) = shedneighbors;
            shedsizes(shedneighbors) = shedsizes(shedneighbors) + 1;
        elseif length(shedneighbors) > 1 % if mulitple, do something else
            % If the point is touching two separated watersheds, decide
            % what to do: if both sheds are large then keep them separated,
            % otherwise if one shed is small and one is large then combine
            if shedsizes(shedneighbors(1)) > size_tol2 && ...
                    shedsizes(shedneighbors(2)) > size_tol2
                shed(ip,jp) = nan;
            else
                shed(ip,jp) = shedneighbors(1);
                shedsizes(shedneighbors(1)) = shedsizes(shedneighbors(1)) ...
                    + 1 + shedsizes(shedneighbors(2));
                shed(shed==shedneighbors(2)) = shedneighbors(1);
                shedsizes(shedneighbors(2)) = 0;
            end
        end
    end
    % Prune off sheds which are too small
    for i=1:shedmax
        if shedsizes(i) < size_tol && shedsizes(i) ~= 0
            shed(shed==i) = nan;
            shedsizes(i) = 0;
        end
    end
    mask = ~isnan(shed);
    FA_results = regionprops(mask, 'PixelList', 'Area', 'Eccentricity', 'Orientation', ...
        'MajorAxisLength', 'MinorAxisLength', 'Centroid');
    
    % Plot the results
    if strcmp(doplot, 'yes')
        figure, imshow(mask, []);
        hold on
        phi = linspace(0,2*pi,50);
        cosphi = cos(phi);
        sinphi = sin(phi);
        for k = 1:length(FA_results)
            xbar = FA_results(k).Centroid(1);
            ybar = FA_results(k).Centroid(2);
            a = FA_results(k).MajorAxisLength/2;
            b = FA_results(k).MinorAxisLength/2;
            theta = pi*FA_results(k).Orientation/180;
            R = [ cos(theta)   sin(theta)
                -sin(theta)   cos(theta)];
            xy = [a*cosphi; b*sinphi];
            xy = R*xy;
            x = xy(1,:) + xbar;
            y = xy(2,:) + ybar;
            plot(x,y,'r','LineWidth',2);
        end
        hold off
    end
else
    FA_results = [];
end
outstruct.FA_results = FA_results;
% End focal adhesion analysis


% OPTIONAL: plot the results
if strcmp(doplot, 'yes')
    % First the weights
    figure('units','pixels',...
      'position',[100 100 size(normim,2)+20 size(normim,1)+20],...
      'name','Weights from Input Image');
    image(S_calc_weights, 'Cdatamapping', 'scaled')
    axis tight
    axis image
    colormap(gray)
    
    % Next the image of orientations
    figure('units','pixels',...
      'position',[100 100 size(orientim,2)+120 size(orientim,1)+20],...
      'name','Orientation Angles');
    pcolor(orientim)
    shading flat
    colorbar
    axis ij
    axis tight
    axis image
    colormap(jet)
    
    
%     % Image of the coherence
%     figure('units','pixels',...
%       'position',[100 100 size(coherence,2)+120 size(coherence,1)+20],...
%       'name','Coherence');
%     pcolor(coherence)
%     shading flat
%     colorbar
%     axis ij
%     axis tight
%     axis image
%     colormap(jet)
    
    
    % Make a color-coded image of nematic order parameter (xy coordinates)
    vals = excise(orientim(:));
    rgbimg = 256 .* (orientim - min(vals)) ./ (max(vals) - min(vals));
    rgbimg = ind2rgb(uint8(rgbimg), jet(256));
    for i=1:size(rgbimg,1)
        for j=1:size(rgbimg,2)
            if isnan(orientim(i,j))
                rgbimg(i,j,1) = .5;
                rgbimg(i,j,2) = .5;
                rgbimg(i,j,3) = .5;
            end
        end
    end
    brightness = 98;
    grayimg = (im - min(im(:))) ./ (prctile(im(:),brightness) - min(im(:)));
    endimg(:,:,1) = rgbimg(:,:,1) .* grayimg;
    endimg(:,:,2) = rgbimg(:,:,2) .* grayimg;
    endimg(:,:,3) = rgbimg(:,:,3) .* grayimg;
    figure('units','pixels',...
      'position',[100 100 size(endimg,2) size(endimg,1)],...
      'resize','off','name','XY Coordinates Orientation Results');
    axes('unit', 'pix', 'position', [1 1 size(endimg,2) size(endimg,1)]);
    image(endimg);
    axis off
    imwrite(endimg,'calculate_order_parameter.bmp');
    
    
    % Make a color-coded image of coherence
%     vals = excise(coherence(:));
%     rgbimg = 256 .* (coherence - min(vals)) ./ (max(vals) - min(vals));
%     rgbimg = ind2rgb(uint8(rgbimg), jet(256));
%     for i=1:size(rgbimg,1)
%         for j=1:size(rgbimg,2)
%             if isnan(coherence(i,j))
%                 rgbimg(i,j,1) = .5;
%                 rgbimg(i,j,2) = .5;
%                 rgbimg(i,j,3) = .5;
%             end
%         end
%     end
%     brightness = 98;
%     grayimg = (im - min(im(:))) ./ (prctile(im(:),brightness) - min(im(:)));
%     endimg(:,:,1) = rgbimg(:,:,1) .* grayimg;
%     endimg(:,:,2) = rgbimg(:,:,2) .* grayimg;
%     endimg(:,:,3) = rgbimg(:,:,3) .* grayimg;
%     figure('units','pixels',...
%       'position',[100 100 size(endimg,2) size(endimg,1)],...
%       'resize','off','name','Coherence Map Results');
%     axes('unit', 'pix', 'position', [1 1 size(endimg,2) size(endimg,1)]);
%     image(endimg);
%     axis off
%     imwrite(endimg,'calculate_order_parameter_coherence.bmp');
    
%     % Make a color-coded image of nematic order parameter (radial
%     % coordinates)
%     vals = excise(radial_orient(:));
%     rgbimg = 2 .* (radial_orient + 1) ./ 2;
%     testmap = [0 .8 0; .7 0 .7];
%     rgbimg = ind2rgb(uint8(rgbimg), testmap);
%     for i=1:size(rgbimg,1)
%         for j=1:size(rgbimg,2)
%             if isnan(radial_orient(i,j))
%                 rgbimg(i,j,1) = .5;
%                 rgbimg(i,j,2) = .5;
%                 rgbimg(i,j,3) = .5;
%             end
%         end
%     end
% %     grayimg = (normim - min(normim(:))) ./ (prctile(normim(:),99) - min(normim(:)));
%     grayimg = (im - min(im(:))) ./ (prctile(im(:),99.6) - min(im(:)));
%     endimg(:,:,1) = rgbimg(:,:,1) .* grayimg;
%     endimg(:,:,2) = rgbimg(:,:,2) .* grayimg;
%     endimg(:,:,3) = rgbimg(:,:,3) .* grayimg;
%     figure('units','pixels',...
%       'position',[100 100 size(endimg,2)+20 size(endimg,1)+20],...
%       'resize','off','name','Radial Coordinates Orientation Results');
%     axes('unit', 'pix', 'position', [10 10 size(endimg,2) size(endimg,1)]);
%     image(endimg);
%     hold on
%     plot(jCOM,iCOM,'rx', 'LineWidth',20)
%     axis off
end
% end data plotting


end
