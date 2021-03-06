% MAIN()                     constructs the Skipping Stone GUI
%
% The GUI is constructed programatically, because the numerous
% checks and hide/show operations and external functions make
% a GUIDE-GUI impractical. Moreover, debugging is easier, there
% are no *.FIG-files necessary (portability), and it is a good
% way to teach new authors the ins-and-outs of how to make a
% GUI.


% Please report bugs and inquiries to:
%
% Name   : Rody P.S. Oldenhuis
% E-mail : oldenhuis@gmail.com
% Licence: 2-clause BSD (See License.txt)


% If you find this work useful, please consider a donation:
% https://www.paypal.me/RodyO/3.5


% If you would like to cite this work, please use the following template:
%
% Rody Oldenhuis, orcid.org/0000-0002-3162-3660. "Skipping Stone" version
% <version>, <date you last used it>. MATLAB interplanetary space mission
% design tool.
% https://nl.mathworks.com/matlabcentral/fileexchange/29272-skipping-stone




% TODO-LIST
%{
     - Use MODEL.AEROGRAVABLE() to check whether the selected body
       can actually be used for an aerograv. If not, generate warning,
       change the selectable bodies to those that CAN be used for
       aerograv and select the first body by default. If yes, just
       change the selectable bodies to aerogravable bodies. Remember
       to change the list back to normal (model.GAMable) when aerograv
       is changed again.

     - Make "General settings" menu, containing "Use MultiThreading",
       "Use NVIDIA-CUDA", and general settings like "Program paths" and
       "File name format" (for BATCH optimizations) etc.

     - "Save as defaults...". this would require an external *.INI file
       for the compiled version (ASCII). This is for someone else
       to do ^_^

     - "Objectives" don't belong on the "Algorithms" tab

     - models should be "clearable" (if model is changed, previous model
        should be unloaded to avoid memory problems)

     - options panes! (for optimizer, exposins, ...)

     - also plot options (number of points, colors, ...)
      (but that's REALLY not important)

 TODO - help and documentation!!
 TODO - CONVERT/Kep2Para: make the program work as in the comments
 TODO - CONVERT/Kep2Para\para2cart - think about efficiency
       (create matrix R = [c au bv] i.s.o. current?)

%}

%% Construct main window, and load all defaults

function varargout = main(varargin) %#ok<VANUS>

    % global variables
    global MainWin
    global launch_tab sequence_tab arrival_tab algorithms_tab output_tab
    global Pareto_tab trajectory_tab central_body_speed post_processing
    global BATCH_optimization optimization_statistics

    % give tab numbers a name for clarity
    % tabs on main window       % tabs on output tab
    launch_tab      = 1;       Pareto_tab         = 1;
    sequence_tab    = 2;       trajectory_tab     = 2;
    arrival_tab     = 3;       central_body_speed = 3;
    algorithms_tab  = 4;       post_processing    = 4;
    output_tab      = 5;       BATCH_optimization = 5;
                               optimization_statistics = 6;
    % Current version
    Version = '0.9 alpha'; % working, but far from perfect...:D

    % return version if requested
    if nargin == 1
        varargout{1} = Version;
        return;
    end

    % default figure handle
    % (figure 123456798 will probably never be used)
    MainWin = 123456789;

    % set the root dir and add GUI subdirectory to MATLAB path
    rootdir = fileparts(mfilename('fullpath'));
    addpath([rootdir, filesep, 'GUI']);

    % load default values
    [environment, model, constants, calculation, settings] = ...
        modify_settings('default_values' , rootdir);

    % Check if another instance is running
    if ishandle(MainWin)
        errordlg({[environment.program_name, ' is already running;']
                   'Only one instance allowed'}, 'One instance is allowed.')
        return
    end

    % Enable multithreaded-BLAS routines
    % NOTE: implicit multithreading only
    % (this is not on by default in MATLAB versions prior to R2009a)
    if verLessThan('MATLAB', '7.8')
        maxNumCompThreads('automatic'); end

    % now build the GUI
    build_main_window(environment,...
                      model,...
                      constants,...
                      calculation,...
                      settings);
end % main
