%% Widefield demo 2
%
% A demo on how to work with widefield data in batch
% (do plab_widefield_demo first)
%
% This code shows how to align widefield data across animals and mice. Note
% that these functions are likely to change/improve over time.
%
% If you would like to use these functions for regular analysis, I would
% suggest that you create a script that has these functions in a pipeline
% to run for your animals.

%% Github repositories to download

% These repositories are necessary for this demo:
%
% https://github.com/PetersNeuroLab/PetersLab_analysis
% https://github.com/petersaj/AP_scripts_peterslab
% https://github.com/kwikteam/npy-matlab

%% Location for saved alignments

% Widefield alignments are saved as a tranform matrix for each animal and
% day. When an alignment is applied (e.g. when loading data), it loads the
% alignment for that animal and applies the transform matrix for that day. 

% Alignments for all animals are kept in: 
% \\qnap-ap001.dpag.ox.ac.uk\APlab\Lab\widefield_alignment
% 
% This has two folders: 
% animal_alignment - the transform matrices for each animal
% retinotopy - retinotopy for each animal used in across-animal aligning
% 
% There are also general files which are not editable, e.g. master_vfs.mat
% is the 'master' retinotopy that all animals are aligned to

% Keeping all alignments in one place makes the pipeline easier, but it
% also means you need to be careful not to mess up other people's data.
% Hopefully this is unlikely, since everyone will only create/use
% alignments for their own animals.


%% Widefield alignment function

% I use the function 'plab.wf.wf_align' to create and apply widefield
% alignments. Check out the help for 'plab.wf.wf_align' for details about how it
% works. Below are examples for how to use it. 

% The general syntax is: 
% im_aligned = wf_align(im_unaligned,animal,day,align_type,master_align)
% 
% This applies alignments to a set of images 'im_unaligned' to create
% aligned images 'im_aligned', or creates and saves new alignments.


%% Creating alignments 

% We can use this animal as an example (data can be overwritten):
animal = 'AP013';

% There are two steps to creating alignments: 

%------- 1) aligning animal across days

% Images are aligned across days using the average blue image. Aligment
% here is rigid, meaning only translation and rotation, since the scaling
% should not change across days. All images are aligned to the first day,
% making the first day the 'master' image for each animal.

% The function automatically finds and loads these, so no unaligned images
% need to be entered (the first argument is blank). This also grabs the
% date for each image, so no dates need to be entered (the third argument
% is blank). Only the second argument (animal), and third argument
% (alignment type) need to be entered. 

% Create a new across-day alignment for this animal by running the command
% below. This will display "Rigid aligning images..." in the command line
% while this works, and then display a scrollable set of images:
plab.wf.wf_align([],animal,[],'new_days');

% The command window should now display "Manual align? (#,s=save, q=quit)".
% entering 's' will save the alignments as-is, entering 'q' will quit
% without saving, and entering a number (e.g. 3) will let you manually
% align that day. 

% Try scrolling through the images, and note that there is a craniotomy on
% day 2 and 3, so they look different from day 1. Days 1-2 should look
% pretty well aligned, but 3 is slightly off. In the command window, type
% '3' and enter.

% This brings up the control point selection tool. Manually click a few
% points on the left (day 3) and right (day 1 = master) images that match
% each other. For the best alignment, select 4+ points, and make the points
% as widely disributed as possible. When you're finished selecting, close
% the window by pressing 'x'. 

% The aligned images will now be displayed again. If you're happy with this
% alignment, you can enter 's'. If you want to do more manual alignment,
% you can enter the number of the day you want to align.

% Once you enter 's', the command window will display "Saved day transforms
% for AP013."


%------- 2) aligning animal to a 'master' alignment (across animals)

% After across-day alignments have been created, you can now align an
% animal to the 'master' alignment, which lets multiple animals be aligned
% to each other. 

% We do this using visual sign maps (vfs), where we align each animal's vfs
% to a 'master' vfs stored in the alignment folder, which was created from
% the average of many aligned mice. Here is what the master vfs looks like:
load('\\qnap-ap001.dpag.ox.ac.uk\APlab\Lab\widefield_alignment\master_vfs.mat');
figure;
imagesc(master_vfs);
axis image;
title('Master VFS');

% In widefield demo 1, we used 'plab.wf.retinotopy_vfs' to create the
% visual field sign map. As a reminder, here we can find sparse noise
% recordings for this animal, and create and plot the vfs. There are more
% than one sparse noise recordings, so this will only plot the first:
recordings = plab.find_recordings(animal,[],'sparse_noise');
plot_recording = 1;
vfs = plab.wf.retinotopy_vfs(animal, ...
    recordings(plot_recording).day,recordings(plot_recording).recording{1});

% To automatically find and save the vfs for all sparse noise recordings
% for a given animal, use the function 'plab.wf.retinotopy_vfs_batch',
% which runs the above retinotopy function and saves into the alignment
% path:
plab.wf.retinotopy_vfs_batch(animal);

% When it's finished, the command window should print the location where it
% saved the vfs maps. 

% After the vfs maps are saved, we can create the across-animal alignment.
% This first creates an average vfs by loading in saved vfs maps, aligning
% across days, and averaging. It then aligns the average vfs for that
% animal with the master vfs. Create that alignment with this command:
plab.wf.wf_align([],animal,[],'new_animal');

% When it's finished running, it will display the master and aligned
% images, and the command window will print "Save animal alignment for
% AP013? (y/n):". If you are happy with the alignment, enter 'y'. If you
% are not, you can exit without saving with 'n'. This doesn't have a manual
% option at the moment, so if this alignment isn't good, then it can't be
% corrected. 

% Enter 'y', and the command window should display "Saved new animal
% alignment for AP013."

% After this step, this animal (AP013) will now have both across-day and
% across-animal alignments.

%% Alignment summary

% In summary, your pipeline for alignment should look like this: 

% Define animal
animal = 'AP013';
% Create across-day alignments
plab.wf.wf_align([],animal,[],'new_days');
% Get and save VFS maps for animal
plab.wf.retinotopy_vfs_batch(animal);
% Create across-animal alignments
plab.wf.wf_align([],animal,[],'new_animal');


%% Applying alignments

% Once an alignment is created, it can be applied with the same function by
% entering the recording day. 

% For example, we can load in an average image from AP013 and align it: 
animal = 'AP013';
rec_day = '2023-11-20';
avg_im_filename = plab.locations.filename('server',animal,rec_day,[],'widefield','meanImage_blue.npy');
avg_im = readNPY(avg_im_filename);
avg_im_aligned = plab.wf.wf_align(avg_im,animal,rec_day);

figure;
subplot(1,2,1);
imagesc(avg_im); axis image;
title('Unaligned image');
subplot(1,2,2);
imagesc(avg_im_aligned); axis image;
title('Aligned image');

% In order to align SVD images, we can align the U's (spatial components)
% in the same way, i.e. 'U_aligned =
% plab.wf.wf_align(U_unaligned,animal,rec_day)'. This is automatically done
% in 'ap.load_experiment' if an alignment exists. 

% Once a widefield image is aligned, you can draw cortical areas on top
% using the function 'ap.wf_draw'. This shows the aligned average image,
% and draws the cortical areas on top: 
figure;
imagesc(avg_im_aligned); axis image;
ap.wf_draw('ccf','r');
title('Aligned image with cortical areas');

%% Using common SVD components across animals

% Using SVD to decompose widefield images is helpful because it it turns
% ~400x400 pixels (160,000 data sources) into ~2000 components, greatly
% reducing the memory requirements when working with the data. 

% Ideally, we would keep the SVD format when analyzing multiple data sets.
% For example, if we wanted to average responses across two days, we would
% rather average the V's (temporal components) than the pixels, because it
% takes up much less space.

% We cannot do this directly, because the U's (spatial components) are
% unique for each recording, so comparing V's across recordings is like
% comparing apples to oranges.

% In order to compare V's, those V's need to be associated with the same
% U's. We can accomplish this with a "change of basis", by creating new V's
% associated with U's of our choice, rather than the original U's. This
% works because there's more than one way to reconstruct the data (or,
% reconstruct it well enough to be useful).

% For example, if you had U = [3,4], trying to reconstruct data point [10],
% you could use V = [2,1]. This gives us U*V = [3,4]*[2,1]' = 10. 
% If we wanted to use a new set of U's U = [9,1], we could reconstruct the
% same data with a new set of V's, e.g. V = [1,1]. This would result in the
% same reconstruction as before, through U*V = [9,1]*[1,1]' = 10.

% In other words, given an SVD output U*V = [data], if we have a set of
% spatial components U_new, we can create a new set of temporal components
% V_new to produce the same reconstructed data U_new*V_new = [data];

% Our lab has a set spatial components 'U_master' which can serve as common
% spatial components across recordings. These are components created from
% the U's of many animals, so are designed to reconstruct data across
% animals. 

% IMPORTANT: these master U's are aligned to the master retinotopy
% described in the alignment section above. This means that you can only
% change the U's from an experiment into the master U's AFTER across-animal
% alignment has been done. Otherwise, the reconstruction won't be very
% good.

% This loads in the master U's and plots them - try scrolling through them
% to see how they might be similar or different to U's you've seen before: 
master_U_fn = fullfile(plab.locations.server_path,'Lab', ...
    'widefield_alignment','U_master.mat');
load(master_U_fn);
ap.imscroll(U_master); 
axis equal;

% Now let's load an example recording:
animal = 'AP019';
rec_day = '2024-04-07';
recodings=plab.find_recordings(animal,rec_day)
rec_time = '0957';
verbose = true;
ap.load_recording;

% This has loaded in 'wf_U' and 'wf_V' as the SVD components from this
% particular recording. 

% [EXERCISE] 
% Look at the master spatial components 'U_master' compared to the
% recording-specific spatial components 'wf_U'. What's the difference? 

% This is the command to convert the SVD components to the U_master: the
% inputs are the recording-specific U and V ('wf_U' and 'wf_V'), and the
% outputs are the master spatial componets 'U_master' and the new
% associated temporal components 'V_master'
% (note: this function loads and outputs 'U_master', so there is no need to
% load it beforehand):
[U_master,V_master] = plab.wf.u2master(wf_U,wf_V);

% [EXERCISE] 
% Reconstruct a collection of frames (e.g. frames 1000-1100) using the
% recording-specific components ('wf_U' and 'wf_V') and the master
% components ('U_master' and 'V_master'). Plot and scroll through them:
% what are the similaries and differences? Do you think the master
% components do a decent job of reconstructing the data produced from the
% recording-specific components?

% [EXERCISE] 
% The loaded workflow is from passive visual stimulus presentations (as in
% the other demos, the stimulus onset times are in 'stimOn_times', and the
% stimulus positions are in 'trial_events.values.TrialStimX'). Create
% averaged stimulus-aligned movies from either the recording-specific
% (wf_V) or master (V_master) components. How well is the data
% reconstructed by the master components?

% [EXERCISE] 
% This animal also had passive visual recordings the next day (2023-08-08).
% Create stimulus-aligned master V's for both the first day (2023-08-07,
% recording 1616) and the second day (2023-08-08, recording 1849). Average
% these together, and make a movie of the average stimulus response across
% these days.


%% Aligning and converting U's on loading

% If you're using my function to load (ap.load_recording, which calls
% ap.load_widefield to load and process the widefield data), this always
% aligns the widefield data, if an alignment exists. 
%
% If you also want to convert the U's to the master Us, you can set a flag
% to do that during loading. 
%
% The variable to set flags during loading is 'load_parts'. These two lines
% set the widefield data to be loaded, and then set the widefield to be
% converted to the master U's;
load_parts = struct;
load_parts.widefield = true;
load_parts.widefield_master = true;

% Now if you load the recording, the 'wf_U' and 'wf_V' variables will be
% converted into the master U's. We can load the data from above, but this
% time it will be aligned/converted during loading: 
animal = 'AP010';
rec_day = '2023-08-07';
rec_time = '1616';
verbose = true;
ap.load_recording;













