function handles = nback_matlab
%% INFORMATION:
% 
% -Keep the 'nback_matlab.m' function within the same folder as the subfolder 
%   'Audio', which contains audio files required for the program to run. As
%   long as 'nback_matlab.m' is on Matlab's path, the program will handle
%   adding the audio files to the path automatically.
% 
% -A 'UserData' subfolder will be written to the same location as
%   'nback_matlab.m the first time the program is run. User performance
%   data will be stored in a .mat file within this folder
%   ('user_data.mat'). Also, whenever the program closes, the
%   'nback_matlab_settings.txt' file will be updated with the user's most
%   recent preferences (this will be loaded upon program reboot).
%   -the user_data.mat file contains a variable 'summaryStats' that
%   summarizes A' (see below), hit rate, false alarm rate, and lure error
%   rate for each stimuli type (position, sound, color); it also contains a
%   'trialData' cell array that shows actual user responses for every trial
%   (rows of table) of each round of n-back played (cells)
% 
% *Scoring:  A' (Discrimination Index), similar to d'
%   H = Hit Rate (hits / # signal trials), 
%   F = False-positive Rate (false pos / # noise trials)
%   aPrime = .5 + sign(H - F) * ((H - F)^2 + abs(H - F)) / (4 * max(H, F) - 4 * H * F);
% 
% A' references:
%   Snodgrass, J. G., & Corwin, J. (1988). Pragmatics of measuring recognition
%       memory: applications to dementia and amnesia. Journal of Experimental
%       Psychology: General, 117(1), 34.
% 
%   Stanislaw, H., & Todorov, N. (1999). Calculation of signal detection theory
%       measures. Behavior research methods, instruments, & computers, 31(1),
%       137-149.
% 
%   *For mathematical formula, see Stanislaw & Todorov (1999, p. 142) 
%  
% *Default Settings:
%   nback = 2; % nBack level
%   percLure = .2; % percentage of non-match trials to be converted to
%       "lure" or interference trials (a lure is n-back +/- 1, requiring
%       engagement of executive control to avoid errors)
%   sound_on = 1; % sound-type nBack: 1 (on), 0 (off)
%   position_on = 1; % position-type nBack: 1 (on), 0 (off)
%   color_on = 0; % color-type nBack: 1 (on), 0 (off)
%   completely_random = 0; % off (0), on (1) (if on, next 3 settings are irrelevant)
%   n_positionHits = 4; % control # of matches for position
%   n_soundHits = 4; % control # of matches for sound
%   n_colorHits = 4;  % control # of matches for color
%   advance_thresh = .9; % minimum score (A') required to advance automatically to next nBack level; 
%   fallback_thresh = .75; % score (A') below which one automatically regresses to previous nBack level; 
%   nTrials = nback + 20; % # trials
%   trial_time = 2.4; % seconds
%   volume_factor = 1; % Default: 1 (no change); >1 = increase volume; <1 = decrease volume
%   sound_type = 1; % use Numbers-Female (1), Numbers-Male (2), Letters-Female1 (3), or Letters-Female2 (4)
%   enable_applause = 1; % 1 (on) or 0 (off); applause on advance to next n-back level
%   enable_boos = 0; % 1 (on) or 0 (off); comical boos if fallback to previous n-back level
%   realtime_feedback = 1; % 1 (on) or 0 (off); highlight labels red/green based on accuracy
%   figure_window_feedback = 1; % 1 (on) or 0 (off); display hits, misses, false alarms in figure window upon completion
%   command_line_feedback = 1; % 1 (on) or 0 (off); display hits, misses, false alarms in command window upon completion
%   show_remaining = 1; % 1 (on) or 0 (off); show remaining trials
%   background_color_scheme = 1; % black (1), white (2)
% 
% *Author: Elliot A. Layden (2017-18), 
% at The Environmental Neuroscience Laboratory, The University of Chicago
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Hotkey Designations: 
startSessionKey = 'space'; 
positionMatchKey = 'a'; 
soundMatchKey = 'l';
colorMatchKey = 'j';
stopRunningKey = 'escape';

%% ------------------------ Begin ------------------------

customizable = {'nback','percLure','sound_on','position_on','color_on','completely_random',...
    'n_soundHits','n_positionHits','n_colorHits','advance_thresh',...
    'fallback_thresh','nTrials','trial_time','volume_factor',...
    'sound_type','enable_applause','enable_boos','realtime_feedback',...
    'command_line_feedback','figure_window_feedback','show_remaining',...
    'background_color_scheme'};

% Initialize Global Variables:
curr_running = 0; stop_running = 0; ix = 1; canResize = false;
position_mem = []; color_mem = []; sound_mem = [];
position_match_vec = []; sound_match_vec = []; color_match_vec = [];
user_position_match = []; user_sound_match = []; user_color_match = [];
position_lures_vec = []; sound_lures_vec = []; color_lures_vec = []; 
advance_thresh = []; fallback_thresh = [];
figure_window_feedback = []; command_line_feedback = []; volume_factor = [];
h_pos_feedback1_pos = []; h_pos_feedback2_pos = [];
h_color_feedback1_pos = []; h_color_feedback2_pos = [];
h_sound_feedback1_pos = []; h_sound_feedback2_pos = [];
percLure = .2; timerVal = 0; times = []; extraKeys = {};

% Initialize object handles:
h_title = 10.384372; h_volume = 28.1717839;
h_pos_feedback1 = 10.384372; h_pos_feedback2 = 10.384372;
h_sound_feedback1 = 10.384372; h_sound_feedback2 = 10.384372;
h_color_feedback1 = 10.384372; h_color_feedback2 = 10.384372;
square_width = .2;
positions = repmat([.05,.05,square_width,square_width],9,1);

% Get script path and subfolders:
script_fullpath = mfilename('fullpath');
[script_path,~,~] = fileparts(script_fullpath);
addpath(genpath(script_path))
audio_path = fullfile(script_path,'Audio');
audio_subfolders = {fullfile(audio_path,'numbers-female'),fullfile(audio_path,'numbers-male'),...
    fullfile(audio_path,'letters-female1'),fullfile(audio_path,'letters-female2')};

% Get Settings:
get_settings;

% Calculate Parameters
[~,nBack_spec] = min(abs((1:20)-nback));
poss_lures = 0:.05:1;
[~,lure_spec] = min(abs(poss_lures-percLure));
trial_num_opts = 15:5:100; % nTrials
n_trial_num_opts = length(trial_num_opts);
if nTrials==nback+20
    trial_num_spec = 1;
else trial_num_spec = 1+min(abs(trial_num_opts-nTrials));
end
trial_time_opts = 1:.2:4; % seconds;
n_trial_time_opts = length(trial_time_opts);
[~,trial_time_spec] = min(abs(trial_time_opts-trial_time));

% Appearance Customization:
if background_color_scheme==1
    txt_color = ones(1,3);
    background_color = zeros(1,3);
    colors = {'y','m','c','r','g','b','w'};
elseif background_color_scheme==2
    txt_color = zeros(1,3);
    background_color = ones(1,3);
    colors = {'y','m','c','r','g','b','k'};
end

% Load Main Audio:
audio_dat = struct('sound',cell(1,4),'Fs',cell(1,4),'nSound',cell(1,4));
for i = 1:4
    listing = dir(audio_subfolders{i});
    nSound = length(listing)-2;
    audio_dat(i).sound = cell(1,nSound);
    audio_dat(i).Fs = cell(1,nSound);
    audio_dat(i).nSound = nSound;
    for j = 1:nSound
        [audio_dat(i).sound{j},audio_dat(i).Fs{j}] = audioread(fullfile(audio_subfolders{i},listing(j+2).name));
    end
end

% Load Reaction Sounds (Applause/Boos):
[applause_dat, Fs_applause] = audioread(fullfile(audio_path,'applause.mp3'));
[boo_dat, Fs_boo] = audioread(fullfile(audio_path,'boo.mp3'));

% Initialize Figure:
handles.figure = figure('menubar','none','color',background_color,'numbertitle',...
    'off','name',['nback_matlab: (',num2str(nback),'-Back)'],'units','norm',...
    'Position',[.271,.109,.488,.802],'ResizeFcn',@resizeScreen); 
handles.axes1 = gca; set(handles.axes1,'Position',[0,0,1,1],'XLim',[0,1],...
    'YLim',[0,1],'Visible','off');
h_title = annotation('textbox','String','','FontName','Helvetica','FontSize',14,...
    'Position',[.4,.955,.2,.05],'Color',txt_color,'HorizontalAlignment','center',...
    'FontWeight','bold','EdgeColor','none');

%% Menus:
% File Menu:
file_menu = uimenu(handles.figure,'Label','File');
%     uimenu(file_menu,'Label','Progress Report','Callback',@show_progress);
    uimenu(file_menu,'Label','Exit','Callback',@close_nback_matlab);
% Session Menu:
session_menu = uimenu(handles.figure,'Label','Session');
    % N-Back Level:
    nBack_num_menu = uimenu(session_menu,'Label','N-Back Level');
    h_nBack = zeros(1,20);
    for i = 1:20
        h_nBack(i) = uimenu(nBack_num_menu,'Label',sprintf('%02g',i),'Callback',{@change_nBack,i});
    end; set(h_nBack(nBack_spec),'Checked','on')
    % Percent Lures:
    percLure_menu = uimenu(session_menu,'Label','Percent Lures/Interference');
    h_lures = zeros(1,21); 
    for i = 1:21
        h_lures(i) = uimenu(percLure_menu,'Label',sprintf('%02g',poss_lures(i)),'Callback',{@change_lures,i});
    end; set(h_lures(lure_spec),'Checked','on')
    % Stimuli Type(s) Menu:
    nBack_type_menu = uimenu(session_menu,'Label','Stimuli Type(s)');
        if position_on
            types_menu(1) = uimenu(nBack_type_menu,'Label','Position','Checked','on','Callback',{@change_types,1});
        else types_menu(1) = uimenu(nBack_type_menu,'Label','Position','Checked','off','Callback',{@change_types,1});
        end
        if sound_on
            types_menu(2) = uimenu(nBack_type_menu,'Label','Sound','Checked','on','Callback',{@change_types,2});
        else types_menu(2) = uimenu(nBack_type_menu,'Label','Sound','Checked','off','Callback',{@change_types,2});
        end
        if color_on
            types_menu(3) = uimenu(nBack_type_menu,'Label','Color','Checked','on','Callback',{@change_types,3});
        else types_menu(3) = uimenu(nBack_type_menu,'Label','Color','Checked','off','Callback',{@change_types,3});
        end
    trial_num_menu = uimenu(session_menu,'Label','# Trials');
        h_trial_num = zeros(1,n_trial_num_opts+1);
        h_trial_num(1) = uimenu(trial_num_menu,'Label','20 + N (default)','Callback',{@change_trial_num,1});
        for i = 2:n_trial_num_opts+1
            h_trial_num(i) = uimenu(trial_num_menu,'Label',sprintf('%1g',trial_num_opts(i-1)),'Callback',{@change_trial_num,i});
        end; set(h_trial_num(trial_num_spec),'Checked','on')
    trial_length_menu = uimenu(session_menu,'Label','Trial Length');
        h_trial_times = zeros(1,n_trial_time_opts);
        for i = 1:n_trial_time_opts
            h_trial_times(i) = uimenu(trial_length_menu,'Label',sprintf('%1.1f seconds',trial_time_opts(i)),'Callback',{@change_trial_time,i});
        end; set(h_trial_times(trial_time_spec),'Checked','on')
advanced_menu = uimenu(session_menu,'Label','Advanced');  
    num_hits_menu = uimenu(advanced_menu,'Label','# Matches');
        sound_matches_menu = uimenu(num_hits_menu,'Label','# Sound Matches');
        n_poss = round(.6*nTrials);
        h_n_sound_matches = zeros(1,n_poss);
        for i = 1:n_poss
            h_n_sound_matches(i) = uimenu(sound_matches_menu,'Label',sprintf('%02g',i),'Callback',{@change_sound_matches,i});
        end
        set(h_n_sound_matches(n_soundHits),'Checked','on')
        position_matches_menu = uimenu(num_hits_menu,'Label','# Position Matches');
        h_n_position_matches = zeros(1,n_poss);
        for i = 1:n_poss
            h_n_position_matches(i) = uimenu(position_matches_menu,'Label',sprintf('%02g',i),'Callback',{@change_position_matches,i});
        end
        set(h_n_position_matches(n_positionHits),'Checked','on')
        color_matches_menu = uimenu(num_hits_menu,'Label','# Color Matches');
        h_n_color_matches = zeros(1,n_poss);
        for i = 1:n_poss
            h_n_color_matches(i) = uimenu(color_matches_menu,'Label',sprintf('%02g',i),'Callback',{@change_color_matches,i});
        end
        set(h_n_color_matches(n_colorHits),'Checked','on')
        if completely_random
            uimenu(num_hits_menu,'Label',...
                'Completely Random','Checked','on','Callback',@change_completely_random);
        else
            uimenu(num_hits_menu,'Label',...
                'Completely Random','Checked','off','Callback',@change_completely_random);
        end
    advance_fallback_menu = uimenu(advanced_menu,'Label','Advance/Fallback');
        advance_menu = uimenu(advance_fallback_menu,'Label','Advance Threshold');
            advance_thresholds = 1:-0.02:0.80; 
            h_advances = zeros(1,length(advance_thresholds));
            for i = 1:length(advance_thresholds)
                h_advances(i) = uimenu(advance_menu,'Label',sprintf('%g%',advance_thresholds(i)),'Callback',{@change_advance_thresh,i});
            end 
            [~,which_advance] = min(abs(advance_thresholds-advance_thresh));
            set(h_advances(which_advance),'Checked','on')
            set(h_advances(6),'Label',[sprintf('%g%',advance_thresholds(6)),' (default)']);
        fallback_menu = uimenu(advance_fallback_menu,'Label','Fallback Threshold');
            fallback_thresholds = .80:-.05:.50; 
            h_fallbacks = zeros(1,length(fallback_thresholds));
            for i = 1:length(fallback_thresholds)
                h_fallbacks(i) = uimenu(fallback_menu,'Label',sprintf('%g%',fallback_thresholds(i)),'Callback',{@change_fallback_thresh,i});
            end
            [~,which_fallback] = min(abs(fallback_thresholds-fallback_thresh));
            set(h_fallbacks(which_fallback),'Checked','on')
            set(h_fallbacks(2),'Label',[sprintf('%g%',fallback_thresholds(2)),' (default)']);
            
feedback_menu = uimenu(handles.figure,'Label','Feedback');
    if show_remaining
        uimenu(feedback_menu,'Label','Show Remaining Trials','Checked','on','Callback',@change_show_remaining);
    else uimenu(feedback_menu,'Label','Show Remaining Trials','Checked','off','Callback',@change_show_remaining);
    end
    if realtime_feedback
        uimenu(feedback_menu,'Label','Realtime Feedback','Checked','on','Callback',@change_realtime_feedback);
    else uimenu(feedback_menu,'Label','Realtime Feedback','Checked','off','Callback',@change_realtime_feedback);
    end
    if figure_window_feedback
        uimenu(feedback_menu,'Label','Results in GUI','Checked','on','Callback',@change_figure_feedback);
    else uimenu(feedback_menu,'Label','Results in GUI','Checked','off','Callback',@change_figure_feedback);
    end
    if command_line_feedback
        uimenu(feedback_menu,'Label','Command Line Results','Checked','on','Callback',@change_cmd_feedback);
    else uimenu(feedback_menu,'Label','Command Line Results','Checked','off','Callback',@change_cmd_feedback);
    end
sound_settings_menu = uimenu(handles.figure,'Label','Sounds');
    uimenu(sound_settings_menu,'Label','Volume','Callback',@change_volume);
    sound_type_menu = uimenu(sound_settings_menu,'Label','Type');
            h_letters = uimenu(sound_type_menu,'Label','Letters');
                h_sound_types(3) = uimenu(h_letters,'Label','Female Voice 1','Callback',{@change_sound_types,3});
                h_sound_types(4) = uimenu(h_letters,'Label','Female Voice 2','Callback',{@change_sound_types,4});
            h_numbers = uimenu(sound_type_menu,'Label','Numbers');  
                h_sound_types(1) = uimenu(h_numbers,'Label','Female Voice','Callback',{@change_sound_types,1});
                h_sound_types(2) = uimenu(h_numbers,'Label','Male Voice','Callback',{@change_sound_types,2});
            set(h_sound_types(sound_type),'Checked','on')
    reactions_menu = uimenu(sound_settings_menu,'Label','Reactions');
        if enable_applause
            uimenu(reactions_menu,'Label','Applause','Checked','on','Callback',@change_applause_setting);
        else uimenu(reactions_menu,'Label','Applause','Checked','off','Callback',@change_applause_setting);
        end
        if enable_boos
            uimenu(reactions_menu,'Label','Boos','Checked','on','Callback',@change_boos_setting);
        else uimenu(reactions_menu,'Label','Boos','Checked','off','Callback',@change_boos_setting);
        end
appearance_menu = uimenu(handles.figure,'Label','Appearance');
    background_menu = uimenu(appearance_menu,'Label','Background');
    if background_color_scheme==1
        h_background_menu1 = uimenu(background_menu,'Label','Black','Checked','on','Callback',{@change_background,1});
        h_background_menu2 = uimenu(background_menu,'Label','White','Callback',{@change_background,2});
        txt_color = ones(1,3);
    elseif background_color_scheme==2
        h_background_menu1 = uimenu(background_menu,'Label','Black','Callback',{@change_background,1});
        h_background_menu2 = uimenu(background_menu,'Label','White','Checked','on','Callback',{@change_background,2});
        txt_color = zeros(1,3);
    else warning('''background_color'' should be either integer 1 (black) or 2 (white)')
        h_background_menu1 = uimenu(background_menu,'Label','Black','Checked','on','Callback',{@change_background,1});
        h_background_menu2 = uimenu(background_menu,'Label','White','Callback',{@change_background,2});
        txt_color = ones(1,3);
    end
   
%% Add other display elements (gridlines, rectangle, etc.):

% Gridlines:
minMeasure = .05; maxMeasure = .95; one_third = .3497; two_thirds = .6503; 
h_grid_lines = zeros(1,4);
h_grid_lines(1) = annotation('line','X',[one_third,one_third],'Y',[minMeasure,maxMeasure],'Color',txt_color,'LineWidth',1);
h_grid_lines(2) = annotation('line','X',[two_thirds,two_thirds],'Y',[minMeasure,maxMeasure],'Color',txt_color,'LineWidth',1);
h_grid_lines(3) = annotation('line','X',[minMeasure,maxMeasure],'Y',[two_thirds,two_thirds],'Color',txt_color,'LineWidth',1);
h_grid_lines(4) = annotation('line','X',[minMeasure,maxMeasure],'Y',[one_third,one_third],'Color',txt_color,'LineWidth',1);

% Initialize Match Button Text and Rectangle:
handles.h_rect = rectangle('Parent',handles.axes1,'Position',...
    [.05,.05,square_width,square_width],'Curvature',[.2,.2],'FaceColor',...
    'blue','Visible','off');
h_txt_begin = annotation('textbox',[.37,.48,.26,.04],'String',...
    ['Press ',startSessionKey,sprintf(' \nto begin')],...
    'Color',txt_color,'HorizontalAlignment','center','VerticalAlignment',...
    'middle','FontSize',14,'FontWeight','bold','EdgeColor','none');
h_txt_pos = text('position',[.04,.01],'String',['Position Match:  ''',positionMatchKey,''''],...
    'Color',txt_color,'HorizontalAlignment','center','VerticalAlignment','top',...
    'FontSize',13,'FontWeight','normal','Visible','off','EdgeColor','none');
h_txt_color = text('position',[.35,.01],'String',['Color Match:  ''',colorMatchKey,''''],...
    'Color',txt_color,'HorizontalAlignment','center','VerticalAlignment','top',...
    'FontSize',13,'FontWeight','normal','Visible','off','EdgeColor','none');
h_txt_sound = text('position',[.658,.01],'String',['Sound Match:  ''',soundMatchKey,''''],...
    'Color',txt_color,'HorizontalAlignment','center','VerticalAlignment','top',...
    'FontSize',13,'FontWeight','normal','Visible','off','EdgeColor','none');
if position_on; set(h_txt_pos,'Visible','on'); end
if color_on;    set(h_txt_color,'Visible','on'); end
if sound_on;    set(h_txt_sound,'Visible','on'); end

% Declare Figure Callbacks:
set(handles.figure,'WindowKeyPressFcn',@keypress_callback);
set(handles.figure,'CloseRequestFcn',@close_nback_matlab)

% Resize Figure function:
canResize = true;
resizeScreen;

%% Run Function:
function run_session
    extraKeys = cell(nTrials,1); curr_running = 1;
    
    % Get Trial Times:
    targetTimes = linspace(1,trial_time * (nTrials-1) + 1, nTrials);
    pause(1 - (toc(timerVal)-0));
    
    disableMenus; % disallow users from toggling settings during play
    if ishandle(h_pos_feedback1)
        delete(h_pos_feedback1); delete(h_pos_feedback2)
    end
    if ishandle(h_sound_feedback1)
        delete(h_sound_feedback1); delete(h_sound_feedback2) 
    end
    if ishandle(h_color_feedback1)
        delete(h_color_feedback1); delete(h_color_feedback2) 
    end
    % START TRIALS:
    ix = 0; times = zeros(nTrials,1);
    while ix<nTrials && ~stop_running
        ix = ix + 1; times(ix) = toc(timerVal);
        
        if realtime_feedback
            set(h_txt_pos,'Color',txt_color)
            set(h_txt_sound,'Color',txt_color)
            set(h_txt_color,'Color',txt_color)
        end
        if show_remaining
            set(h_title,'String',num2str(nTrials-ix))
        end
        tic
        if position_on
            set(handles.h_rect,'Position',positions(position_mem(ix),:),'Visible','on')
        end
        if color_on
            set(handles.h_rect,'FaceColor',colors{color_mem(ix)})
        end
        if sound_on
            sound(audio_dat(sound_type).sound{sound_mem(ix)}*volume_factor,audio_dat(sound_type).Fs{sound_mem(ix)})
        end
        
        if abs(toc(timerVal)-targetTimes(ix)) < .5
            pause(trial_time - (toc(timerVal)-targetTimes(ix)));
        else
            pause(trial_time - sign(toc(timerVal)-targetTimes(ix))*.5);
        end
    end
    curr_running = 0; set(handles.h_rect,'Visible','off')
    set(h_txt_begin,'Visible','on'); set(h_title,'String','')
    set(h_txt_pos,'Color',txt_color); set(h_txt_sound,'Color',txt_color)
    set(h_txt_color,'Color',txt_color)
    % Calculate Score:
    if ~stop_running
        
        % Position Stats:
        position_hits = sum((position_match_vec+user_position_match)==2);
        position_false = sum((position_match_vec-user_position_match)==-1);
        position_lure_errors = sum((position_lures_vec+user_position_match)==2)/sum(position_lures_vec);
        H_pos = position_hits / n_positionHits; F_pos = position_false / (nTrials - n_positionHits);
        pos_score = .5 + sign(H_pos - F_pos) * ((H_pos - F_pos)^2 + abs(H_pos - F_pos)) / (4 * max(H_pos, F_pos) - 4 * H_pos * F_pos);
        
        position_stats = cell(5,2); 
        position_stats(1:5,1) = {'Position','A'':','Hits:','False-Alarms:','Lure Errors:'};
        position_stats{2,2} = num2str(round(pos_score,2));
        position_stats{3,2} = [num2str(round(100*H_pos)),'%'];
        position_stats{4,2} = [num2str(round(100*F_pos)),'%'];
        position_stats{5,2} = [num2str(round(100*position_lure_errors)),'%'];
        
        % Sound Stats:
        sound_hits = sum((sound_match_vec+user_sound_match)==2);
        sound_false = sum((sound_match_vec-user_sound_match)==-1);
        sound_lure_errors = sum((sound_lures_vec+user_sound_match)==2)/sum(sound_lures_vec);
        H_sound = sound_hits / n_soundHits; F_sound = sound_false / (nTrials - n_soundHits);
        sound_score = .5 + sign(H_sound - F_sound) * ((H_sound - F_sound)^2 + abs(H_sound - F_sound)) / (4 * max(H_sound, F_sound) - 4 * H_sound * F_sound);
        
        sound_stats = cell(5,2); 
        sound_stats(1:5,1) = {'Sound','A'':','Hits:','False-Alarms:','Lure Errors:'};
        sound_stats{2,2} = num2str(round(sound_score,2));
        sound_stats{3,2} = [num2str(round(100*H_sound)),'%'];
        sound_stats{4,2} = [num2str(round(100*F_sound)),'%'];
        sound_stats{5,2} = [num2str(round(100*sound_lure_errors)),'%'];
        
        % Color Stats:
        color_hits = sum((color_match_vec+user_color_match)==2);
        color_false = sum((color_match_vec-user_color_match)==-1);
        color_lure_errors = sum((color_lures_vec+user_color_match)==2)/sum(color_lures_vec);
        H_color = color_hits / n_colorHits; F_color = color_false / (nTrials - n_colorHits);
        color_score = .5 + sign(H_color - F_color) * ((H_color - F_color)^2 + abs(H_color - F_color)) / (4 * max(H_color, F_color) - 4 * H_color * F_color);
        
        color_stats = cell(5,2); 
        color_stats(1:5,1) = {'Color','A'':','Hits:','False-Alarms:','Lure Errors:'};
        color_stats{2,2} = num2str(round(color_score,2));
        color_stats{3,2} = [num2str(round(100*H_color,1)),'%'];
        color_stats{4,2} = [num2str(round(100*F_color,1)),'%'];
        color_stats{5,2} = [num2str(round(100*color_lure_errors)),'%'];
     
        % Command-line Feedback:
        if command_line_feedback
            if position_on && sound_on && color_on 
                session_stats = [position_stats,sound_stats,color_stats];
            elseif position_on && sound_on && ~color_on 
                session_stats = [position_stats,sound_stats];
            elseif position_on && ~sound_on && ~color_on
                session_stats = position_stats;
            elseif ~position_on && sound_on && ~color_on 
                session_stats = sound_stats;
            elseif ~position_on && sound_on && color_on 
                session_stats = [sound_stats,color_stat];
            elseif position_on && ~sound_on && color_on  
                session_stats = [position_stats,color_stats];
            elseif ~position_on && ~sound_on && color_on 
                session_stats = color_stats;
            else return;
            end
            disp(session_stats)
        end
        
        % Write data:
        user_position_match = user_position_match==1;
        user_sound_match = user_sound_match==1;
        user_color_match = user_color_match==1;
        writeData(pos_score, position_lure_errors, H_pos, F_pos, ...
            sound_score, sound_lure_errors, H_sound, F_sound, ...
            color_score, color_lure_errors, H_color, F_color, ...
            (1:nTrials)', times, position_match_vec', position_lures_vec', ...
            user_position_match', position_mem', sound_match_vec', sound_lures_vec', ...
            user_sound_match', sound_mem', color_match_vec', color_lures_vec', ...
            user_color_match', color_mem',extraKeys)
        
        % Check if Advance:
        advance_vec = [(pos_score >= advance_thresh),(sound_score >= advance_thresh),(color_score >= advance_thresh)];
        fallback_vec = [(pos_score <= fallback_thresh),(sound_score <= fallback_thresh),(color_score <= fallback_thresh)];
        if all(advance_vec(logical([position_on,sound_on,color_on])))
            if command_line_feedback
                disp('Congratulations, you have been advanced to the next nBack level!')
            end
            nback = nback+1;
            if get(h_trial_num(1),'checked')
                nTrials = nTrials + 1;
            end
            set(handles.figure,'name',['nback_matlab: (',num2str(nback),'-Back)'])
            set(h_nBack(nback-1),'Checked','off'); set(h_nBack(nback),'Checked','on');
            if enable_applause; sound(applause_dat,Fs_applause); end
            if figure_window_feedback
                figure_feedback_color = 'g';
            end
        elseif any(fallback_vec(logical([position_on,sound_on,color_on])))
            if nback>1; 
                nback = nback-1; 
                if get(h_trial_num(1),'checked')
                    nTrials = nTrials - 1;
                end
                if command_line_feedback
                    disp('Good effort - nBack level lowered.')
                end
                set(handles.figure,'name',['nback_matlab: (',num2str(nback),'-Back)'])
                set(h_nBack(nback+1),'Checked','off'); set(h_nBack(nback),'Checked','on');
            else if command_line_feedback; disp('Good effort.'); end
            end
            if enable_boos; sound(boo_dat,Fs_boo); end
            if figure_window_feedback
                figure_feedback_color = 'r';
            end
        else
            if command_line_feedback; disp('Good score - keep training!'); end
            if figure_window_feedback
                if background_color_scheme==1
                    figure_feedback_color = 'w';
                elseif background_color_scheme==2
                    figure_feedback_color = 'k';
                end
            end
        end
        if figure_window_feedback
            if position_on
                main_char = char('Position','A'':','Hits:','False-Alarms:','Lure Errors:');
                char2 = char('',position_stats{2,2},position_stats{3,2},position_stats{4,2},position_stats{5,2});
                h_pos_feedback1 = text('String',main_char,'FontName','Helvetica','FontSize',13,...
                    'Position',h_pos_feedback1_pos,'Color',figure_feedback_color,...
                    'FontWeight','bold','EdgeColor','none','HorizontalAlignment','left');
                h_pos_feedback2 = text('String',char2,'FontName','Helvetica','FontSize',13,...
                    'Position',h_pos_feedback2_pos,'Color',figure_feedback_color,'HorizontalAlignment','right',...
                    'FontWeight','normal','EdgeColor','none');
            end
            if sound_on
                main_char = char('Sound','A'':','Hits:','False-Alarms:','Lure Errors:');
                char2 = char('',sound_stats{2,2},sound_stats{3,2},sound_stats{4,2},sound_stats{5,2});             
                h_sound_feedback1 = text('String',main_char,'FontName','Helvetica','FontSize',13,...
                    'Position',h_sound_feedback1_pos,'Color',figure_feedback_color,...
                    'FontWeight','bold','EdgeColor','none');
                h_sound_feedback2 = text('String',char2,'FontName','Helvetica','FontSize',13,...
                    'Position',h_sound_feedback2_pos,'Color',figure_feedback_color,'HorizontalAlignment','right',...
                    'FontWeight','normal','EdgeColor','none');
            end
            if color_on
                main_char = char('Color','A'':','Hits:','False-Alarms:','Lure Errors:');
                char2 = char('',color_stats{2,2},color_stats{3,2},color_stats{4,2},color_stats{5,2});
                h_color_feedback1 = text('String',main_char,'FontName','Helvetica','FontSize',13,...
                    'Position',h_color_feedback1_pos,'Color',figure_feedback_color,...
                    'FontWeight','bold','EdgeColor','none');
                h_color_feedback2 = text('String',char2,'FontName','Helvetica','FontSize',13,...
                    'Position',h_color_feedback2_pos,'Color',figure_feedback_color,'HorizontalAlignment','right',...
                    'FontWeight','normal','EdgeColor','none'); 
            end
        end
    else stop_running = 0;
    end
    enableMenus;
end

% Key Press Callback:
function keypress_callback(~,which_key,~)
    switch which_key.Key
        case startSessionKey
            if ~curr_running
                timerVal = tic;
                curr_running = 1; set(h_txt_begin,'Visible','off')
                if completely_random
                    position_mem = randi(9,[1,nTrials]);
                    color_mem = randi(7,[1,nTrials]);
                    sound_mem = randi(audio_dat(sound_type).nSound,[1,nTrials]);
                else
                    % Generate Idx:
                    wasSuccess = false;
                    while ~wasSuccess
                        [position_mem, position_match_vec, position_lures_vec, wasSuccess] = generateIdx(nTrials, nback, n_positionHits, percLure, 9);
                    end; wasSuccess = false;
                    while ~wasSuccess
                        [sound_mem, sound_match_vec, sound_lures_vec, wasSuccess] = generateIdx(nTrials, nback, n_soundHits, percLure, audio_dat(sound_type).nSound);
                    end; wasSuccess = false;
                    while ~wasSuccess
                        [color_mem, color_match_vec, color_lures_vec, wasSuccess] = generateIdx(nTrials, nback, n_colorHits, percLure, 7);
                    end
                    user_position_match = nan(1, nTrials); 
                    user_sound_match = nan(1, nTrials); 
                    user_color_match = nan(1, nTrials); 
                end
                run_session
            end
        case positionMatchKey
            if curr_running
                user_position_match(ix) = 1;
                if realtime_feedback
                    if position_match_vec(ix)
                        set(h_txt_pos,'Color',[0,1,0])
                    else
                        set(h_txt_pos,'Color',[1,0,0]) 
                    end
                end
            end
        case colorMatchKey
            if curr_running
                user_color_match(ix) = 1;
                if realtime_feedback
                    if color_match_vec(ix)
                        set(h_txt_color,'Color',[0,1,0])
                    else
                        set(h_txt_color,'Color',[1,0,0]) 
                    end
                end
            end
        case soundMatchKey
            if curr_running
                user_sound_match(ix) = 1;
                if realtime_feedback
                    if sound_match_vec(ix)
                        set(h_txt_sound,'Color',[0,1,0])
                    else
                        set(h_txt_sound,'Color',[1,0,0]) 
                    end
                end
            end
        case stopRunningKey
            stop_running = 1;
        otherwise
            if curr_running
                extraKeys{ix} = [extraKeys{ix}, ',', which_key.Key];
            end
    end
end

% Change nBack Types:
function change_types(~,~,which_type)
    switch which_type
        case 1,
            switch position_on
                case 1, set(types_menu(1),'Checked','off'); position_on = false;
                    set(h_txt_pos,'Visible','off')
                case 0, set(types_menu(1),'Checked','on'); position_on = true;
                    set(h_txt_pos,'Visible','on')
            end
        case 2,
            switch sound_on
                case 1, set(types_menu(2),'Checked','off'); sound_on = false;
                    set(h_txt_sound,'Visible','off')
                case 0, set(types_menu(2),'Checked','on'); sound_on = true;
                    set(h_txt_sound,'Visible','on')
            end
        case 3,
            switch color_on
                case 1, set(types_menu(3),'Checked','off'); color_on = false;
                    set(h_txt_color,'Visible','off')
                    set(handles.h_rect,'FaceColor','b')
                case 0, set(types_menu(3),'Checked','on'); color_on = true;
                    set(h_txt_color,'Visible','on')
            end
    end                 
end
    
% Change nBack Level:
function change_nBack(~,~,nBack_spec)
    nback = nBack_spec;
    for ix1 = 1:20
        set(h_nBack(ix1),'Checked','off');
    end
    set(h_nBack(nback),'Checked','on')
    set(handles.figure,'name',['nback_matlab: (',num2str(nback),'-Back)'])
end

function change_lures(~,~,lure_spec)
    percLure = poss_lures(lure_spec);
    for ix1 = 1:21
        set(h_lures(ix1),'Checked','off');
    end
    set(h_lures(lure_spec),'Checked','on')
end

% Change N Position Matches:
function change_position_matches(~,~,change_n)
    n_positionHits = change_n;
    for ixx = 1:n_poss
        set(h_n_position_matches(ixx),'Checked','off');
    end
    set(h_n_position_matches(n_positionHits),'Checked','on')
end

% Change N Sound Matches:
function change_sound_matches(~,~,change_n)
    n_soundHits = change_n;
    for ixx = 1:n_poss
        set(h_n_sound_matches(ixx),'Checked','off');
    end
    set(h_n_sound_matches(n_soundHits),'Checked','on')
end

% Change N Color Matches:
function change_color_matches(~,~,change_n)
    n_colorHits = change_n;
    for ixx = 1:n_poss
        set(h_n_color_matches(ixx),'Checked','off');
    end
    set(h_n_color_matches(n_colorHits),'Checked','on')
end

% Change Completely Random:
function change_completely_random(hObject,~,~)
    switch completely_random
        case 1, set(hObject,'Checked','off'); completely_random = 0;
        case 0, set(hObject,'Checked','on'); completely_random = 1;
    end
end

% Change Advance Threshold:
function change_advance_thresh(~,~,advance_spec)
	advance_thresh = advance_thresholds(advance_spec);
    for ixx = 1:length(advance_thresholds)
        set(h_advances(ixx),'Checked','off')
    end
    set(h_advances(advance_spec),'Checked','on')
end

% Change Advance Threshold:
function change_fallback_thresh(~,~,fallback_spec)
	fallback_thresh = fallback_thresholds(fallback_spec);
    for ixx = 1:length(fallback_thresholds)
        set(h_fallbacks(ixx),'Checked','off')
    end
    set(h_fallbacks(fallback_spec),'Checked','on')
end

% Change # Trials:
function change_trial_num(~,~,trial_num_spec)
    if (trial_num_spec > 1)
        nTrials = trial_num_opts(trial_num_spec-1);
    else 
        nTrials = nback + 20;
    end
    for ixx = 1:n_trial_num_opts+1
        set(h_trial_num(ixx),'Checked','off')
    end
    set(h_trial_num(trial_num_spec),'Checked','on')
end

% Change Trial Time:
function change_trial_time(~,~,trial_time_spec)
	trial_time = trial_time_opts(trial_time_spec);
    for ixx = 1:n_trial_time_opts
        set(h_trial_times(ixx),'Checked','off')
    end
    set(h_trial_times(trial_time_spec),'Checked','on')
end

% Change Show Remaining Trials:
function change_show_remaining(hObject,~,~)
    switch show_remaining
        case 1, set(hObject,'Checked','off'); show_remaining = 0;
        case 0, set(hObject,'Checked','on'); show_remaining = 1;
    end
end

function change_realtime_feedback(hObject,~,~)
    switch realtime_feedback
        case 1, set(hObject,'Checked','off'); realtime_feedback = 0;
        case 0, set(hObject,'Checked','on'); realtime_feedback = 1;
    end
end

function change_figure_feedback(hObject,~,~)
    switch figure_window_feedback
        case 1, set(hObject,'Checked','off'); figure_window_feedback = 0;
        case 0, set(hObject,'Checked','on'); figure_window_feedback = 1;
    end
end

function change_cmd_feedback(hObject,~,~)
    switch command_line_feedback
        case 1, set(hObject,'Checked','off'); command_line_feedback = 0;
        case 0, set(hObject,'Checked','on'); command_line_feedback = 1;
    end
end

% Change Sound Settings:
function change_volume(~,~,~)
    if ~ishandle(h_volume)
    h_volume = figure('menubar','none','color',background_color,'numbertitle',...
        'off','name','Adjust Volume','units','norm','Position',[.4729,.8646,.2174,.0352]);
    uicontrol(h_volume,'Style','slider','units','normalized','Position',...
        [0,.2,1,.6],'Min',0,'Max',1.2,'Value',volume_factor,'Callback',@change_volume_slider);
    else figure(h_volume)
    end
end

function change_volume_slider(hObject,~,~)
    volume_factor = hObject.Value;
end

function change_sound_types(~,~,which_type)
    for ixx = 1:4; set(h_sound_types(ixx),'Checked','off'); end
    set(h_sound_types(which_type),'Checked','on')
    sound_type = which_type;
end

function change_applause_setting(hObject,~,~)
    switch enable_applause
        case 1, set(hObject,'Checked','off'); enable_applause = 0;
        case 0, set(hObject,'Checked','on'); enable_applause = 1;
    end
end

function change_boos_setting(hObject,~,~)
    switch enable_boos
        case 1, set(hObject,'Checked','off'); enable_boos = 0;
        case 0, set(hObject,'Checked','on'); enable_boos = 1;
    end
end

% Appearance:
function change_background(~,~,which_background)
    switch which_background
        case 1, 
            if strcmp(h_background_menu1.Checked,'on')
                background_color_scheme = 2;
                colors = {'y','m','c','r','g','b','k'};
                set(h_background_menu1,'Checked','off')
                set(h_background_menu2,'Checked','on')
            elseif strcmp(h_background_menu1.Checked,'off')
                background_color_scheme = 1;
                colors = {'y','m','c','r','g','b','w'};
                set(h_background_menu1,'Checked','on')
                set(h_background_menu2,'Checked','off')
            end
        case 2,
            if strcmp(h_background_menu2.Checked,'on')
                background_color_scheme = 1;
                colors = {'y','m','c','r','g','b','w'};
                set(h_background_menu1,'Checked','on')
                set(h_background_menu2,'Checked','off')
            elseif strcmp(h_background_menu2.Checked,'off')
                background_color_scheme = 2;
                colors = {'y','m','c','r','g','b','k'};
                set(h_background_menu1,'Checked','off')
                set(h_background_menu2,'Checked','on')
            end
    end           
    if background_color_scheme==1
        txt_color = ones(1,3);
        set(handles.figure,'Color',zeros(1,3))
    elseif background_color_scheme==2
        txt_color = zeros(1,3);
        set(handles.figure,'Color',ones(1,3))
    end
    for ixxx = 1:4
        set(h_grid_lines(ixxx),'Color',txt_color)
    end
    set(h_title,'Color',txt_color); set(h_txt_begin,'Color',txt_color)
    set(h_txt_pos,'Color',txt_color); set(h_txt_color,'Color',txt_color)
    set(h_txt_sound,'Color',txt_color)
end

% Write Current Trial Data:
function writeData(pos_score, pos_lure_errors, H_pos, F_pos, ...
        sound_score, sound_lure_errors, H_sound, F_sound, ...
        color_score, color_lure_errors, H_color, F_color, trials, times, ...
        position_matches, pos_lures, position_user, position_stimuli,...
        sound_matches, sound_lures, sound_user, sound_stimuli, ...
        color_matches, color_lures, color_user, color_stimuli, extraKeys)
    
        if ~isdir(fullfile(script_path,'UserData')); mkdir(fullfile(script_path,'UserData')); end
        dataPath = fullfile(script_path,'UserData','user_data.mat');
        
        if exist(dataPath,'file')
            load(dataPath);
            summaryStats2 = table({char(datetime)}, nback, pos_score, ...
                pos_lure_errors, H_pos, F_pos, sound_score, ...
                sound_lure_errors, H_sound, F_sound, color_score, ...
                color_lure_errors, H_color, F_color,'VariableNames',...
                {'date_time','nback_level',...
                'A_pos','lure_errors_pos','hits_pos','false_alarms_pos',...
                'A_sound','lure_errors_sound','hits_sound','false_alarms_sound',...
                'A_color','lure_errors_color','hits_color','false_alarms_color'});
            summaryStats = [summaryStats; summaryStats2]; %#ok
            trialData{length(trialData)+1} = table(trials, times, ...
                position_matches,pos_lures, position_user, position_stimuli, ...
                sound_matches, sound_lures, sound_user, sound_stimuli, ...
                color_matches, color_lures, color_user, color_stimuli, extraKeys,...
                'VariableNames',{'trial_num','seconds_from_start',...
                'position_matches','position_lures','user_position_matches','position_stimuli',...
                'sound_matches','sound_lures','user_sound_matches','sound_stimuli',...
                'color_matches','color_lures','user_color_matches','color_stimuli','extra_keys'}); %#ok        
        else
            summaryStats = table({char(datetime)}, nback, pos_score, ...
                pos_lure_errors, H_pos, F_pos, sound_score, ...
                sound_lure_errors, H_sound, F_sound, color_score, ...
                color_lure_errors, H_color, F_color,'VariableNames',...
                {'date_time','nback_level',...
                'A_pos','lure_errors_pos','hits_pos','false_alarms_pos',...
                'A_sound','lure_errors_sound','hits_sound','false_alarms_sound',...
                'A_color','lure_errors_color','hits_color','false_alarms_color'}); %#ok
            trialData{1} = table(trials, times, ...
                position_matches,pos_lures, position_user, position_stimuli, ...
                sound_matches, sound_lures, sound_user, sound_stimuli, ...
                color_matches, color_lures, color_user, color_stimuli, extraKeys,...
                'VariableNames',{'trial_num','seconds_from_start',...
                'position_matches','position_lures','user_position_matches','position_stimuli',...
                'sound_matches','sound_lures','user_sound_matches','sound_stimuli',...
                'color_matches','color_lures','user_color_matches','color_stimuli','extra_keys'}); %#ok       
        end
        save(dataPath,'summaryStats','trialData')
end

% Get Default Settings:
function get_defaults
    nback = 2; % nBack level
    percLure = .2; % proportion of non-match trials to become lures
    nTrials = nback + 20;
    trial_time = 2.4; % seconds
    sound_on = true; % sound-type nBack: 1 (on), 0 (off)
    position_on = true; % position-type nBack: 1 (on), 0 (off)
    color_on = false; % color-type nBack: 1 (on), 0 (off)
    completely_random = 0; % off (0), on (1) (if on, next 3 settings are irrelevant)
    n_positionHits = 4; % control # of matches for position
    n_soundHits = 4; % control # of matches for sound
    n_colorHits = 4;  % control # of matches for color
    advance_thresh = .9; % A' default
    fallback_thresh = .75; % A' default
    volume_factor = 1; % Default: 1 (no change); >1 = increase volume; <1 = decrease volume
    sound_type = 3; % use Numbers-Female (1), Numbers-Male (2), Letters-Female1 (3), or Letters-Female2 (4)
    enable_applause = 1; % 1 (on) or 0 (off)
    enable_boos = 0; % 1 (on) or 0 (off); comical boos if fail to advance to next nBack
    realtime_feedback = 1; % 1 (on) or 0 (off); highlight labels red/green based on accuracy
    figure_window_feedback = 1; % 1 (on) or 0 (off); display hits, misses, false alarms in figure window upon completion
    command_line_feedback = 1; % 1 (on) or 0 (off); display hits, misses, false alarms in command window upon completion
    show_remaining = 1; % 1 (on) or 0 (off); show remaining trials as title
    background_color_scheme = 1; % black (1), white (2)
end

% Get Custom Settings:
function success = get_settings
    if ~isdir(fullfile(script_path,'UserData')); mkdir(fullfile(script_path,'UserData')); end
    listing = dir(fullfile(script_path,'UserData','nback_matlab_settings.txt'));
    if isempty(listing); 
        success = 0;
        get_defaults
        return;
    end
    nback_matlab_settings = importdata(fullfile(script_path,'UserData',listing(1).name));
    for ix1 = 1:length(customizable); eval(nback_matlab_settings{ix1}); end
    success = 1;
end

% Write Custom Settings:
function write_settings
    % Change Settings Data:
    nback_matlab_settings = cell(length(customizable),1);
    for ix1 = 1:length(customizable)
        setting1 = eval(customizable{ix1});
        nback_matlab_settings{ix1} = [customizable{ix1},'=',num2str(setting1),';'];  
    end
    nback_matlab_settings = char(nback_matlab_settings);
    % Write Text File:
    if ~isdir(fullfile(script_path,'UserData')); mkdir(fullfile(script_path,'UserData')); end
    fileID = fopen(fullfile(script_path,'UserData','nback_matlab_settings.txt'),'w');
    for ix1 = 1:length(customizable)
        fprintf(fileID,'%s\r\n',nback_matlab_settings(ix1,:));
    end
    fclose(fileID);
end

% Close Requests:
function close_nback_matlab(varargin)
    % Save Custom Settings:
    write_settings
    % Close:
    delete(handles.figure)
end

function enableMenus
    set(handles.figure,'resize','on')
    set(file_menu,'Enable','on'); set(session_menu,'Enable','on')
    set(nBack_type_menu,'Enable','on'); set(feedback_menu,'Enable','on')
    set(sound_settings_menu,'Enable','on'); set(appearance_menu,'Enable','on') 
end

function disableMenus
    clear sound % stop any sound if playing
    set(handles.figure,'resize','off')
    set(file_menu,'Enable','off')
    set(session_menu,'Enable','off')
    set(nBack_type_menu,'Enable','off')
    set(feedback_menu,'Enable','off')
    set(sound_settings_menu,'Enable','off')
    set(appearance_menu,'Enable','off')
end

function resizeScreen(varargin)    
    if canResize 
        % Determine which dimension needs shrinking to become square:
        for ixx = 1:4; set(h_grid_lines(ixx),'units','norm'); end
        set(h_grid_lines(1),'Y',[minMeasure,maxMeasure])
        set(h_grid_lines(2),'Y',[minMeasure,maxMeasure])
        set(h_grid_lines(3),'X',[minMeasure,maxMeasure])
        set(h_grid_lines(4),'X',[minMeasure,maxMeasure])

        % Change units to points and shrink width or height to match:
        for ixx = 1:4; set(h_grid_lines(ixx),'units','points'); end
        ySpec = get(h_grid_lines(1),'Y'); xSpec = get(h_grid_lines(3),'X');
        height = diff(ySpec); width = diff(xSpec);
        axesRatio = height/width;
        diffPixels = abs(height - width);
        if height > width % shrink height to match
            set(h_grid_lines(1),'Y',[ySpec(1) + diffPixels/2, ySpec(2) - diffPixels/2])
            set(h_grid_lines(2),'Y',[ySpec(1) + diffPixels/2, ySpec(2) - diffPixels/2])
        elseif height < width % shrink width to match
            set(h_grid_lines(3),'X',[xSpec(1) + diffPixels/2, xSpec(2) - diffPixels/2])
            set(h_grid_lines(4),'X',[xSpec(1) + diffPixels/2, xSpec(2) - diffPixels/2])
        end
        
        % Change units back to normalized, and adjust grid spacing:
        for ixx = 1:4; set(h_grid_lines(ixx),'units','norm'); end
        xSpec = get(h_grid_lines(3),'X'); ySpec = get(h_grid_lines(1),'Y'); 
        oneThirdX = (diff(xSpec)/3) + xSpec(1);
        twoThirdsX = xSpec(2) - (diff(xSpec)/3);
        oneThirdY = (diff(ySpec)/3) + ySpec(1);
        twoThirdsY = ySpec(2) - (diff(ySpec)/3);
        set(h_grid_lines(1),'X',[oneThirdX, oneThirdX])
        set(h_grid_lines(2),'X',[twoThirdsX, twoThirdsX])
        set(h_grid_lines(3),'Y',[oneThirdY, oneThirdY])
        set(h_grid_lines(4),'Y',[twoThirdsY, twoThirdsY])

        % Square Positions:
        ylim(handles.axes1,[0,axesRatio]); % adjust Y-axes limits based on its ratio w/ X
        square_width = .95*(twoThirdsX-oneThirdX);
        useGapX = ((twoThirdsX-oneThirdX)-square_width)/2;
        oneThirdY = oneThirdY*axesRatio; twoThirdsY = twoThirdsY*axesRatio;
        useGapY = ((twoThirdsY-oneThirdY)-square_width)/2;
        positions = repmat([.05,.05,square_width,square_width],9,1);
        positions(1,1:2) = [xSpec(1)+useGapX,twoThirdsY+useGapY];
        positions(2,1:2) = [oneThirdX+useGapX,twoThirdsY+useGapY];
        positions(3,1:2) = [twoThirdsX+useGapX,twoThirdsY+useGapY];
        positions(4,1:2) = [xSpec(1)+useGapX,oneThirdY+useGapY];
        positions(5,1:2) = [oneThirdX+useGapX,oneThirdY+useGapY];
        positions(6,1:2) = [twoThirdsX+useGapX,oneThirdY+useGapY];
        positions(7,1:2) = [xSpec(1)+useGapX,ySpec(1)*axesRatio+useGapY];
        positions(8,1:2) = [oneThirdX+useGapX,ySpec(1)*axesRatio+useGapY];
        positions(9,1:2) = [twoThirdsX+useGapX,ySpec(1)*axesRatio+useGapY];
        
        % Adjust match button text positions:
        set(h_txt_pos,'position',[median([xSpec(1),oneThirdX]),ySpec(1)*axesRatio])
        set(h_txt_color,'position',[median([oneThirdX, twoThirdsX]),ySpec(1)*axesRatio])
        set(h_txt_sound,'position',[median([twoThirdsX,xSpec(2)]),ySpec(1)*axesRatio])
        
        % Adjust Feedback Positions:
        h_pos_feedback1_pos = [xSpec(1)+.005*xSpec(1), median([ySpec(1)*axesRatio, oneThirdY])];
        h_pos_feedback2_pos = [oneThirdX-.005*oneThirdX, median([ySpec(1)*axesRatio, oneThirdY])];
        h_color_feedback1_pos = [oneThirdX+.005*oneThirdX, median([ySpec(1)*axesRatio, oneThirdY])];
        h_color_feedback2_pos = [twoThirdsX-.005*twoThirdsX, median([ySpec(1)*axesRatio, oneThirdY])];
        h_sound_feedback1_pos = [twoThirdsX+.005*twoThirdsX, median([ySpec(1)*axesRatio, oneThirdY])];
        h_sound_feedback2_pos = [xSpec(2)-.005*xSpec(2), median([ySpec(1)*axesRatio, oneThirdY])];
        if ishandle(h_pos_feedback1); set(h_pos_feedback1,'position',h_pos_feedback1_pos); end
        if ishandle(h_pos_feedback2); set(h_pos_feedback2,'position',h_pos_feedback2_pos); end
        if ishandle(h_color_feedback1); set(h_color_feedback1,'position',h_color_feedback1_pos); end
        if ishandle(h_color_feedback2); set(h_color_feedback2,'position',h_color_feedback2_pos); end
        if ishandle(h_sound_feedback1); set(h_sound_feedback1,'position',h_sound_feedback1_pos); end
        if ishandle(h_sound_feedback2); set(h_sound_feedback2,'position',h_sound_feedback2_pos); end
    end
end

function [idx, matches, lures, countMatches, countLures, wasSuccess] = generateIdx(nTrials, nback, nMatches, percLure, nStimuli)

    wasSuccess = true;
    
    nPassLimit = 30;
    if nback==1
        nLures = round(percLure * (nTrials - nMatches));
    else
        nLures = round(percLure * (nTrials - nMatches - nback + 1));
    end

    %% 1. Add matches:
    idx = zeros(1,nTrials); nMatchesReal = 0; nLuresReal = 0; nPasses = 0;
    while (nMatchesReal ~= nMatches || nLuresReal > nLures) && nPasses<nPassLimit
        nPasses = nPasses + 1;
        matchInd = find(idx==0); matchInd(matchInd<=nback) = [];
        matchInd = matchInd(randperm(length(matchInd),length(matchInd)));
        if nMatchesReal < nMatches
            if  idx(matchInd(1)-nback)>0 % match already present
                idx(matchInd(1)) = idx(matchInd(1)-nback);
            else % match not present
                idx(matchInd(1)) = randi(nStimuli);
                idx(matchInd(1)-nback) = idx(matchInd(1));
            end
        elseif ((nMatchesReal > nMatches) && (nLuresReal <= nLures))
            opts = find(idx>0); 
            if ~isempty(opts)
                idx(opts(randi(length(opts),1))) = 0;
                idx = clearNonMatches(idx, nback);
            end
        elseif (nLuresReal > nLures)
            [~, lureInd] = findLures(idx,nback);
            idx(lureInd) = 0;
        end
        nMatchesReal = findMatches(idx, nback); nLuresReal = findLures(idx, nback);
    end
    idx = clearNonMatches(idx, nback);
    
    %% 2. Add Lures:
    lurePoss = find(idx==0);
    if nback>1; 
        lurePoss(lurePoss<nback) = [];
    else lurePoss(lurePoss<nback+1) = [];
    end
    lurePossShuffled = lurePoss(randperm(length(lurePoss),length(lurePoss)));

    [~, whereMatches] = findMatches(idx, nback);
    whereMatches = find(whereMatches);
    whereMatches = [whereMatches, whereMatches-nback]; 
    
    nPasses = 0; ixx = 0; nLuresReal = findLures(idx, nback);
    while nLuresReal~=nLures && ixx<length(lurePossShuffled) && nPasses < nPassLimit
        nPasses = nPasses + 1;
        if nLuresReal < nLures
            % Iterate:
            ixx = ixx + 1;
            lureInd = lurePossShuffled(ixx);
            
            % Check if viable index:
            if idx(lureInd)~=0; continue; end
            
            % Choose Lure Type (1 or 2):
            onlyOneOption = false;
            if lureInd-nback-1 <= 0 || any(whereMatches == lureInd-nback-1)
                lureType = 2; onlyOneOption = true;
                if lureInd-nback+1<=0
                    continue; % can't add this lure at all 
                end
            elseif lureInd-nback+1 <= 0 || any(whereMatches == lureInd-nback+1)
                continue;
            else
                lureType = randi(2);
            end
            
            % Find stimuli to avoid at this index:
            toAvoid = whichToAvoid(lureInd, lureType, idx, nback);
            indOpts = setdiff(1:nStimuli, toAvoid);
            useStimulus = indOpts(randi(length(indOpts),1));
            
            if (lureType==1) 
                if (idx(lureInd-nback-1)==0) 
                    idx(lureInd) = useStimulus;
                    idx(lureInd-nback-1) = useStimulus;
                elseif (idx(lureInd-nback-1)~=idx(lureInd-nback)) && (lureInd+nback>nTrials || idx(lureInd-nback-1)~=idx(lureInd+nback)) % if -nBack-1 already has a lure, just use that stimulus for new lure:
                    idx(lureInd) = idx(lureInd-nback-1);
                elseif onlyOneOption
                    continue;
                else  % try adding -nBack+1 instead
                    toAvoid = whichToAvoid(lureInd, 2, idx, nback); 
                    indOpts = setdiff(1:nStimuli, toAvoid);
                    useStimulus = indOpts(randi(length(indOpts),1));
                    if (idx(lureInd-nback+1)==0) 
                        idx(lureInd) = useStimulus;
                        idx(lureInd - nback + 1) = useStimulus;
                    elseif ((lureInd-nback<=0 || idx(lureInd - nback + 1) ~= idx(lureInd - nback)) && (lureInd+nback>nTrials || idx(lureInd - nback + 1) ~= idx(lureInd + nback))) 
                            idx(lureInd) = idx(lureInd - nback + 1);
                    else 
                        continue;
                    end
                end
            elseif (lureType==2) 
                if (idx(lureInd-nback+1)==0) 
                    idx(lureInd) = useStimulus;
                    idx(lureInd - nback + 1) = useStimulus;
                elseif ((lureInd-nback<=0 || idx(lureInd - nback + 1) ~= idx(lureInd - nback)) && (lureInd+nback>nTrials || idx(lureInd - nback + 1) ~= idx(lureInd + nback)))
                    idx(lureInd) = idx(lureInd - nback + 1);
                elseif (onlyOneOption) 
                    continue;
                else  % Try adding -nBack-1 instead:
                    toAvoid = whichToAvoid(lureInd, 1, idx, nback); 
                    indOpts = setdiff(1:nStimuli, toAvoid);
                    useStimulus = indOpts(randi(length(indOpts),1));
                    if (idx(lureInd-nback-1)==0) 
                        idx(lureInd) = useStimulus;
                        idx(lureInd - nback - 1) = useStimulus;
                    elseif ((lureInd-nback-1<=0 || idx(lureInd - nback - 1) ~= idx(lureInd-nback)) && (lureInd+nback>nTrials || idx(lureInd - nback - 1) ~= idx(lureInd + nback))) 
                        idx(lureInd) = idx(lureInd - nback - 1);
                    else 
                        continue;
                    end
                end
            end
            
            % Count up the # of lures (new lures can be created unexpectedly):
            nLuresReal = findLures(idx, nback);
        
        elseif nLuresReal > nLures
            opts = lurePossShuffled(1:ixx);
            if ~isempty(opts)
                idx(opts(randi(length(opts),1))) = 0;
            end
            
            % Count up the # of lures (new lures can be created unexpectedly):
            nLuresReal = findLures(idx, nback);
        end
    end
    
    % If still too few lures, add more wherever gaps remain:
    remainingPossible = find(idx==0); nPasses = 0;
    while nLuresReal~=nLures && ~isempty(remainingPossible) && nPasses<nPassLimit 
        nPasses = nPasses + 1;
        if nLuresReal < nLures
            for ixx = remainingPossible
                possible = [];
                if (ixx+nback-1 <= nTrials) 
                    if ((ixx+nback>nTrials || idx(ixx+nback-1)~=idx(ixx+nback)) && (ixx-nback<=0 || idx(ixx+nback-1)~=idx(ixx-nback))) 
                        possible = [possible, idx(ixx+nback-1)]; %#ok
                    end
                    if (ixx+nback+1 <= nTrials) 
                        if (idx(ixx+nback+1)~=idx(ixx+nback) && (ixx-nback<=0 || idx(ixx+nback+1)~=idx(ixx-nback))) 
                            possible = [possible, idx(ixx+nback+1)]; %#ok
                        end
                    end
                end
                if (ixx-nback+1 > 0) 
                    if ((ixx-nback<=0 || idx(ixx-nback+1)~=idx(ixx-nback)) && (ixx+nback>nTrials || idx(ixx-nback+1)~=idx(ixx+nback))) 
                        possible = [possible, idx(ixx-nback+1)]; %#ok
                    end
                    if (ixx-nback-1 > 0) 
                        if (idx(ixx-nback-1)~=idx(ixx-nback) && (ixx+nback>nTrials || idx(ixx-nback-1)~=idx(ixx+nback))) 
                            possible = [possible, idx(ixx-nback-1)]; %#ok
                        end
                    end
                end
                possible(possible==0) = []; %#ok
                if ~isempty(possible)
                    idx(ixx) = possible(randi(length(possible),1));
                    nLuresReal = findLures(idx, nback);
                    remainingPossible = find(idx==0);
                    if (nLuresReal >= nLures) 
                        break;
                    end
                else % if is empty, cannot add lure
                    remainingPossible(remainingPossible==ixx) = [];
                end
            end
        elseif nLuresReal > nLures
            opts = lurePossShuffled(1:min(ixx,length(lurePossShuffled)));
            if ~isempty(opts)
                idx(opts(randi(length(opts),1))) = 0;
            end
            
            % Count up the # of lures (new lures can be created unexpectedly):
            nLuresReal = findLures(idx, nback);
            remainingPossible = find(idx==0);
        end
    end
    
    %% 3. Add non-lures & non-matches:
    for ixx = find(idx==0)
        toAvoid = [];
        if (ixx-nback+1 > 0) 
            if (idx(ixx-nback+1) ~= 0) 
                toAvoid = [toAvoid, idx(ixx-nback+1)]; %#ok
            end
            if (ixx-nback > 0) 
                if (idx(ixx-nback) ~= 0) 
                    toAvoid = [toAvoid, idx(ixx-nback)]; %#ok
                end
                if (ixx-nback-1 > 0) 
                    if (idx(ixx-nback-1) ~= 0) 
                        toAvoid = [toAvoid, idx(ixx-nback-1)]; %#ok
                    end
                end
            end
        end
        if (ixx+nback-1 <= nTrials) 
            if (idx(ixx+nback-1) ~= 0) 
                toAvoid = [toAvoid, idx(ixx+nback-1)]; %#ok
            end
            if (ixx+nback <= nTrials) 
                if (idx(ixx+nback) ~= 0) 
                    toAvoid = [toAvoid, idx(ixx+nback)]; %#ok
                end
                if (ixx+nback+1 <= nTrials) 
                    if (idx(ixx+nback+1) ~= 0) 
                        toAvoid = [toAvoid, idx(ixx+nback+1)]; %#ok
                    end
                end
            end
        end
        indOpts = setdiff(1:nStimuli, toAvoid);
        idx(ixx) = indOpts(randi(length(indOpts),1));
    end
    
    %% Check final counts of matches and lures:
    [countMatches, matches] = findMatches(idx, nback);
    if countMatches~=nMatches; 
        warning(['# matches, target : ',num2str(countMatches),', ',num2str(nMatches)]); 
        wasSuccess = false;
    end
    [countLures, lures] = findLures(idx, nback);
    if countLures~=nLures && ~isempty(remainingPossible) %&& nPasses < nPassLimit
        warning(['# lures, target : ',num2str(countLures),', ',num2str(nLures)]); 
        wasSuccess = false;
        return; 
    end
    
    %% Helper Functions:
    function [nMatches, matches] = findMatches(vec, nback)
        matches = false(size(vec)); nMatches = 0;
        for ixxx = 1:length(vec)
           if vec(ixxx)>0 && ixxx>nback && vec(ixxx)==vec(ixxx-nback)
               matches(ixxx) = true;
               nMatches = nMatches + 1;
           end
        end
    end

    function [nLures, lures] = findLures(vec,nback)
        nLures = 0; lures = false(size(vec));
        for ixxx = 1:length(vec)
            if vec(ixxx)~=0 
                if (ixxx-nback<=0 || vec(ixxx)~=vec(ixxx-nback)) % assure not a match
                    if (ixxx-nback+1>0 && vec(ixxx)==vec(ixxx-nback+1))
                        nLures = nLures+1;
                        lures(ixxx) = true;
                        continue;
                    end
                    if ixxx-nback-1>0 && vec(ixxx)==vec(ixxx-nback-1)
                        nLures = nLures+1;
                        lures(ixxx) = true;
                        continue;
                    end
                end
            end
        end
    end

    function avoid = whichToAvoid(lureInd, whichLureType, idx, nBackNum)
        avoid = [];
        if (lureInd+nBackNum <= length(idx)) && idx(lureInd+nBackNum)~=0 % Avoid new match forward
            avoid = [avoid, idx(lureInd+nBackNum)];
        end
        if (lureInd-nBackNum>0) && idx(lureInd-nBackNum)~=0  % Avoid new match backward
            avoid = [avoid, idx(lureInd-nBackNum)];
        end
        if (whichLureType==1) 
            if (lureInd-nBackNum-1-nBackNum>0) && idx(lureInd-nBackNum-1-nBackNum)~=0  % Avoid creating new match
                avoid = [avoid, idx(lureInd-nBackNum-1-nBackNum)];
            end
            if (idx(lureInd-nBackNum-1+nBackNum) ~= 0) % Avoid creating new match
                avoid = [avoid, idx(lureInd-nBackNum-1+nBackNum)]; 
            end
        elseif (whichLureType==2) 
            if (lureInd-nBackNum-nBackNum+1 > 0) && idx(lureInd-nBackNum-nBackNum+1)~=0  % Avoid creating new match
                avoid = [avoid, idx(lureInd-nBackNum-nBackNum+1)];
            end
            if (lureInd+1<=length(idx)) && idx(lureInd+1)~=0  % Avoid creating new match
                avoid = [avoid, idx(lureInd+1)];
            end
        end
    end

    function vec = clearNonMatches(vec, nback)
        for ixxx = 1:length(vec)
           if vec(ixxx)>0 && ~((ixxx>nback && vec(ixxx)==vec(ixxx-nback)) || (ixxx+nback<=length(vec) && vec(ixxx)==(vec(ixxx+nback))))
               vec(ixxx) = 0;
           end
        end
    end
end

end