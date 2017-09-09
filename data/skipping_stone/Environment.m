classdef Environment
    
    properties (Constant)
      

        program_name = 'SkippingStone';

        % AUTHORS & CONTRIBUTORS
        % more authors can easily be added this way
        % just copy-paste this below, and fill in your own info:
        authors = {'Rody P.S. Oldenhuis',...                      % name
                   'oldenhuis@gmail.com',...                      % contact info
                   'Delft University of Technology, Holland,',... % affiliation field 1
                   'Partly performed at ISAS/JAXA, Japan.',...    % affiliation field 2
                   ' '};                                          % affiliation field 3

        contributors = {'John d`Errico',...
                        'Yi Cao',...
                        'Stephen Morris',...
                        'Erik Johnson',...
                        'ioxv.4623 (all MATLAB FEX)',...
                        'Dr. Dario Izzo (ESA/ACT)',...
                        };

        % Using Octave or MATLAB?        
        ui = octave_or_matlab();

        % pathing: set rootdir, datadir, previous path, and set proper path
        %addpath(genpath(rootdir));
        pathing = struct('rootdir'    , '',...
                         'datadir'    , fullfile(rootdir, 'data'),...
                         'prevpath'   , path(),...
                         'MP_filename', fullfile(ep.datadir, 'asteroids', 'MPCORB.DAT') ...
                         );

        
        colors = struct('window_bgcolor', get(0, 'defaultUicontrolBackgroundColor'),... % background color for windows        
                        'edit_bgcolor'  , 'White'); % Default BG-color for edit boxes should be WHITE

        % Check if the optimization toolbox is available
        optim_toolbox_available = ~isempty(ver('optim'));
       
        % Check for ephemerides
        have_DE405_MEX = exist('JPL_DE405_ephemerides','file')==3;       
        have_DE421_MEX = exist('JPL_DE421_ephemerides','file')==3;

    end
    
end

% Determine whether Octave or MATLAB is running
% (trick from author 'ioxv.4623' on matlab file exchange)
function ui = octave_or_matlab()
    
    ui  = 'Unknown';
    lic = license('inuse');
    
    if any(strcmpi({lic.feature}, 'Matlab'))
        ui = 'Matlab'; end
    if any(strcmpi({lic.feature}, 'Octave'))
        ui = 'Octave'; end
    
end
