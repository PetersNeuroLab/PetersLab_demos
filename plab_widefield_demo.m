%% Widefield demo
%
% A demo on how to work with widefield data
% (do plab_data_demo first)
% Whenever you see a line of code, run it and continue reading
% [EXERCISE] marks a section with instructions to complete

%% Github repositories to download

% These repositories are necessary for this demo:
% https://github.com/PetersNeuroLab/PetersLab_analysis
% https://github.com/petersaj/AP_scripts_peterslab

%% Example dataset

% My loading script (ap.load_recording) loads and preprocesses widefield
% data (specifically in ap.load_widefield), so that will be used for this
% demo.

% Load this data for the following demo: 
animal = 'AP010';
rec_day = '2023-08-07';
rec_time = '1616';
verbose = true;
ap.load_recording;

%% Widefield data format

% Widefield data is captured in pixels at ~400x400, which is 160k pixels
% per frame. This can be hard to work with, because it takes up a lot of
% memory and operations across this many pixels can be slow or impossible.
% For example, the loaded data is just 5 minutes of imaging, and the raw
% widefield images are 8.6GB. The raw values are also 16-bit numbers, but
% many matlab operations require 32-bit, so that requires 8.6GB * 2 =
% 17.2GB RAM just for this short dataset.

% To make the data more accessible, we compress it with singular value
% decomposition (SVD). Instead of each pixel being independent, SVD
% expresses the data as weighted maps of all pixels (U, spatial components)
% which vary over time (V spatial components). The data can then be
% reconstructed by multiplying corresponding U's and V's and summing across
% components, which can be done with matrix multiplication data = U*V.

% The first component explains the most variance in the data, and
% subsequent components explain progressively less variance. This means we
% can either use all components (= number of pixels or timepoints,
% whichever's smaller) to exactly reconstruct the data, or a subset of
% components to reconstruct most of the variance. We keep 2000 components,
% since this explains >95% of the variance. This brings the memory in our
% example from 8.6GB to 1.5GB.

% ap.load_recording loads these components as: 
% wf_U: U (spatial components - Y pixels x X pixels x N components)
% wf_V: V (temporal components - N components x M timepoints)
% wf_t: timestamps in seconds for each frame (M timepoints x 1, in Timelite clock)

% Take a look at the U's (spatial components). Note the early components
% have structure (patterns of pixels that are correlated/anti-correlated
% and explain much of the variance - the first component is the average by
% definition), and later components look like noise (these pick up small
% variance):
ap.imscroll(wf_U);
clim([-0.01,0.01]); 
axis image
colormap(gray);

% Each page of U is a component, which corresponds to the row in V. For
% example, this plots the U and V for component 3:
plot_component = 3;

figure; colormap(gray)
tiledlayout(1,2,'tilespacing','tight');

nexttile;
imagesc(wf_U(:,:,plot_component));
axis image
title(sprintf('U component %d',plot_component))
nexttile; 
plot(wf_t,wf_V(plot_component,:));
title(sprintf('V component %d',plot_component))

% The U's (spatial) and V's (temporal) are discussed above, but SVD
% produces a third matrix S (amplitude of each component), such that
% standard SVD expresses data D as D = USV. For simplicity, we save V's as
% S*V. Because of this, our early V components have larger amplitudes than
% late components because they drive more variance. Here's an example of
% amplitude differences between component 3 and 500, note that 3 is larger
% than 500:
plot_components = [3,500];
figure;
plot(wf_t,wf_V(plot_components,:));
xlabel('Time (s)');
title(sprintf('V components: %d,%d',plot_components(1),plot_components(2)));

% Pixel values for each frame can be reconstructed by matrix multiplication
% of U's and V's (this requires some reshapes from how we normally store
% them: we store U as Y x X pixels x components, but for matrix
% multiplcation it has to be flattened into pixels x components). Here's an
% example of reconstructing one frame:
use_frame = 500;
wf_U_flat = reshape(wf_U,[],size(wf_U,3)); % flatten U for U*V operation
example_fluorescence_flat = wf_U_flat*wf_V(:,use_frame); % matrix multiply all components for given frame
example_fluorescence_frame = reshape(example_fluorescence_flat,size(wf_U,[1,2])); % reshape flat pixels into frame size
figure;
imagesc(example_fluorescence_frame);
axis image off;
title('Example frame fluorescence');

% Note that the image above looks different from the raw image, because it
% is mean-subtracted. The average raw image is also saved and loaded, and
% can be displayed like this:
figure;
imagesc(wf_avg);
axis image off;
title('Average widefield image');

% As mentioned above: the benefit of SVD is that you don't need to use all
% components, because you can capture most of the variance with a subset.
% For example, here we can reconstruct the same frame above, but with
% different numbers of components. Note that more components give the image
% more focused regions of activity, but there isn't much difference in the
% images after you have ~50 components:

% (set the numbers of components to use for each reconstruction)
use_components = [5,10,50,100,200,500];
figure;
tiledlayout(1,length(use_components));
for curr_reconstruct = 1:length(use_components)
    % (set the number of components for this reconstruction)
    curr_components = use_components(curr_reconstruct);
    % (this is as above, but indexing the columns of wf_U_flat - the rows
    % are pixel, and the rows of wf_V - the columns are frame)
    example_fluorescence_flat = wf_U_flat(:,1:curr_components)*wf_V(1:curr_components,use_frame);
    example_fluorescence_frame = reshape(example_fluorescence_flat,size(wf_U,[1,2]));

    nexttile;
    imagesc(example_fluorescence_frame);
    clim([-0.006,0.006])
    axis image off
    title(sprintf('%d components',curr_components));
end

% To avoid writing out the matrix multiplication each time, we have the
% function 'plab.wf.svd2px', which does the reshape and multiplication:
use_frame = 500;
example_fluorescence_frame = plab.wf.svd2px(wf_U,wf_V(:,use_frame));
figure;
imagesc(example_fluorescence_frame);
axis image off;
title('Example frame fluorescence');

% Multiple frames can be reconstructed together (or all of them, if you
% don't index the V's). Scroll through this data to see a movie of the
% reconstructed widefield frames:
use_frames = 300:400;
example_fluorescence_frames = plab.wf.svd2px(wf_U,wf_V(:,use_frames));
ap.imscroll(example_fluorescence_frames);
axis image;

% Another benefit of SVD is that we can do any linear operations (+,-,/,*)
% on the V temporal components, rather than reconstructed pixels, because
% SVD reconstruction is a linear operation. For example, if we wanted to
% average frames together, we can average the V's first and then
% reconstruct, which is equivalent to reconstructing pixels and then
% averaging: 
frames_to_average = 1:10:100; % select frames to average
avg_V = mean(wf_V(:,frames_to_average),2); % average V selected frames for all components
avg_px = plab.wf.svd2px(wf_U,avg_V);
figure;imagesc(avg_px)
axis image off

% [EXERCISE] 
% 1) Show that the above point is true: average frames by (1) averaging
% V's, then reconstructing to (2) reconstructing to pixels, then averaging.
% Write code to check if 1 and 2 are equivalent.
%
% 2) Show that the order of operations does matter for non-linear
% operations: do the same as above but with standard devation (which is
% non- linear because it uses a square root)
%
% 3) The Bonsai workflow in the loaded data is 'lcr_passive', which is
% passive stimulus presentations as described in the demo 'plab_data_demo'.
% Using what you learned from that demo, make a stimulus-triggered average
% movie of widefield fluorescence for -10:30 frames around stimulus
% presentations. Do this separately for stimulus X positions = -90 (left),
% 0 (center) and 90 (right).

% In the same way that we can do linear operations in time on the V's
% above, we can do linear operations in space on the U's. For example, if
% we want to average pixels together in a region-of-interest (ROI), we can
% do that first on the U's and then reconstruct a single trace. 

% Use this code to draw an ROI on the average image:
figure; colormap(gray);
imagesc(wf_avg); 
axis image off
roi_poly = drawpolygon;
% Then use this code to turn your ROI polygon into a binary mask: 
roi_mask = createMask(roi_poly);
figure; imagesc(roi_mask); axis image;
title('ROI mask');

% [EXERCISE] 
% 1) Using the mask 'roi_mask' created above, average the pixel values
% within in the ROI for each component in U, giving you one number for each
% component (rather than a spatial component Y pixels x X pixels x N
% components, you'll have 1 pixel-average x N components). Reconstruct the
% pixel trace with this new averaged spatial component in the same way the
% full-image reconstructions were done above.
%
% 2) Show that this the above method is equivalent to first reconstructing,
% then averaging pixels.

% I have a function 'ap.wf_roi' to make ROIs to make that process easier,
% here's an example use, which can output both the trace and the mask:
[roi_trace,roi_mask2] = ap.wf_roi(wf_U,wf_V,wf_avg);
figure;plot(roi_trace);
title('ROI fluorescence');

% You can also enter a mask as a 5th argument, for example this will use
% exactly the mask you created above, and the trace should be the same as
% the one you reconstructed above:
roi_trace = ap.wf_roi(wf_U,wf_V,wf_avg,[],roi_mask);
figure;plot(roi_trace);
title('ROI fluorescence (with previous mask)');

% This is also an interesting tool to play with: this gives a map of
% correlations between a clicked pixel and all other pixels. Try clicking
% around the brain to get a sense for the pattern of correlations, also you
% can press 'h' to enable 'hover mode' that follows your cursor: 
ap.wf_corrviewer(wf_U,wf_V);

% [EXERCISE] 
% In the correlation viewer above: 
%
% 1) There should be some regions that have topographic relationships,
% where adjacent pixels in one region are correlated to adjacent pixels in
% another region. You should see the hotspots of correlation move in an
% ordered way across multiple regions, and then come together at specific
% points. Which regions have this relationship? Where do the correlations
% join together? 
% 
% 2) If you had to divide the cortex into regions based on the correlation
% patterns, where would you draw the lines? Take a look at Figure 3 in this
% paper which labels dorsal cortical areas:
% https://www.cell.com/cell/fulltext/S0092-8674(20)30402-5.
% What regions do your divisions correspond to?








