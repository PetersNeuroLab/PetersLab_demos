%% Data demo
%
% A demo on how to use the server and load data
% Whenever you see a line of code, run it and continue reading
% [EXERCISE] marks a section with instructions to complete

%% Repos to download

% Clone these repos: 
% https://github.com/PetersNeuroLab/PetersLab_analysis
% https://github.com/petersaj/AP_scripts_peterslab

%% File structure and finding recordings

% To connect to the server: 
% - Open up the file explorer, select This PC in the left hand menu 
% - Right click on this PC and select map network drive 
% - Enter the following address: \\qnap-ap001.dpag.ox.ac.uk\APlab\ 
% - It should prompt you for credentials, please enter your MSD details in the following format: 
% - Username – MSD\(MSD username) 
% - Password – normal MSD password 

% The server contains these top folders:
% \Data: for all raw and preprocessed data (all users can write)
% \Users: all users have a folder for their own use (given users can write)
% \Lab: files for the general lab and rigs (only Andy can write)

% Data is organized on the server in these folders: 
% \Data
% |- animal name (initials & number, e.g. AP012)
% | |- recording date (yyyy-mm-dd, e.g. 2023-03-16)
% |  |- recording time (Recording_HHMM, e.g. Recording_1021 started at 10:21am)
% |  | |- recording data by modality (e.g. mousecam, widefield)
% |  |- day data spanning multiple recordings (e.g. ephys, widefield)
% |- non-recording data (e.g. histology)


%% Finding data from recordings

% The function 'plab.find_recordings' is used to find data from recordings
% The syntax is: plab.find_recordings(animal,recording_day,workflow)
% (type 'help plab.find_recordings' for documentation).
% 'animal' must be filled in, other specifics can be left out to return all
% relevant recordings.

% --- Working with plab.find_recordings

% This line finds all recordings done with animal AP010:
% (day and workflow are not entered)
recordings = plab.find_recordings('AP010');

% 'recordings' is a structure with one entry for each day containing:
% - day: day of recording
% - index: index of recording within the day for that animal (e.g. 3rd recording = [3])
% - recording: time of recording (HHMM)
% - workflow: Bonsai workflow (i.e. stimuli or task protocol)
% - mousecam: if there was mouse camera recording (0/1 = no/yes)
% - widefield: if there was widefield recording (0/1 = no/yes)
% - ephys: if there was electrophysiology recording (0/1 = no/yes)

% [EXERCISE]
% 1) How many days were recordings performed on AP010?
% 2) How many recordings were performed on AP010 across all days?
% 3) How many unique Bonsai workflows were used, and what were they?

% This line finds all recordings of AP010 on 2023-08-24: 
% (workflow is not entered)
recordings = plab.find_recordings('AP010','2023-08-24');

% [EXERCISE]
% Which recording modalities were recorded this day?

% This line finds all recordings of AP010 with workflow 'lcr_passive':
% (date is not entered)
recordings = plab.find_recordings('AP010',[],'lcr_passive');

% [EXERCISE]
% 1) How many recordings of this workflow included widefield? 
% 2) How many recordings of this workflow included electrophysiology? 
% 3) Make a new variable which is a subset of 'recordings' with widefield

% This line finds recordings of AP010 with workflow 'lcr_passive' on
% 2023-08-16 (all fields are entered)
recordings = plab.find_recordings('AP010','2023-08-16','lcr_passive');

% [EXERCISE]
% This day has 2 recordings of the same workflow. This usually is
% because there was an issue with the first one and it needed to be re-run.
% In the number order of that day's recordings, which ones were
% 'lcr_passive'? (hint: 'index')

% The 'workflow' can include * as a wildcard. For example, in the task
% workflow 'stim_wheel_right', there is a 'stage1' and 'stage2' version.
% This returns specifically recordings with stage1:
recordings = plab.find_recordings('AP010',[],'stim_wheel_right_stage1');
% And this returns recordings with either stage: 
recordings = plab.find_recordings('AP010',[],'stim_wheel_right*');

% [EXERCISE]
% Write code to return the day AP010 switched from stage1 to stage2.


%% Constructing server paths

% The class 'plab.locations' is used to find or construct paths on the server.

% --- Standardized generic locations

% This is the general path to the server:
plab.locations.server_path
% This is the data path on the server:
plab.locations.server_data_path
% This is where local data should be kept on each computer: 
plab.locations.local_data_path

% [EXERCISE]
% Use the function 'fullfile' and the above information to generate the
% path: server/Users/(your name)/test_path

% --- Recording locations

% The method 'plab.locations.filename' is used to construct paths/filenames
% for recording data
% (type 'help plab.locations.filename' for documentation). 
% The syntax is:   constructed_filename = filename('server' or 'local',animal,rec_day,rec_time,folder1,...,folderN)

% For example: 
% This constructs the path for animal AP001 on day 2000-01-01 at 12:00:
plab.locations.filename('server','AP001','2000-01-01','1200')
% Input arguments can be added or removed as necessary, for example:
% This constructs the path for the above recording and the subfolder
% 'mousecam'
plab.locations.filename('server','AP001','2000-01-01','1200','mousecam')
% This constructs the path for the above animal and day in the folder
% 'ephys' (not in a recording time folder):
plab.locations.filename('server','AP001','2000-01-01',[],'ephys')
% And this example constructs a path 'histology' in the animal folder (not
% in a day folder): 
plab.locations.filename('server','AP001',[],[],'histology')

% Note that paths with plab.locations.filename are constructed whether or
% not the folder/file exists (e.g. above example paths do not exist).

% [EXERCISE] 
% Use 'plab.find_recording' to find the recording time for AP010 on
% 2023-08-10 with workflow 'lcr_passive', then use
% 'plab.locations.filename' to construct the path to the 'widefield'
% subfolder within the folder for that recording. Write a line to check
% whether that path exists on the server.

%% Loading data

% There isn't standardized lab code for loading data, since this is often
% customized for each person depending on their needs. This block demos my
% code, which can be used as a template, or as-is if it works for you (note
% that it's subject to regular changes).

% My loading code is in my repository (petersaj/AP_scripts_peterslab) and
% is called with 'ap.load_recording' after defining animal/day/time. This
% loads an example dataset with passive stimuli, widefield, and ephys:
animal = 'AP005';
rec_day = '2023-06-21';
rec_time = '1859';
verbose = true; % this turns on/off progress display in command line
ap.load_recording;

% --- Recording locations





% TO LOAD PARTS
% load_parts.widefield = true;
% load_parts.ephys = true;
% load_parts.mousecam = true;



















