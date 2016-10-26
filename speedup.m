% SPEEDUP       Compile several functions to speed up Skipping Stone
%
% Compared to MATLAB 7 (and earlier), later versions of MATLAB execute
% M-code significantly faster thanks to the introduction of the
% JIT-compiler. Thanks to JIT technology, many of the speed-related
% problems associated with older versions of MATLAB have now completely
% or partially disappeared. M-code written in C++/FORTRAN-style (a very
% bad idea before the introduction of JIT) can now also run without
% too much overhead. For many programs, this translates into M-code
% runtimes that are only ~5% slower than compiled C++/FORTRAN-code.
%
% However, JIT-technology does not *always* work; it has a quite
% constrictive set of rules that the code should obey before JIT-compilation
% can actually be applied. Also, at the current stage of development, JIT
% is not very intelligent, and will in many cases incur substantial a
% overhead, reducing the performance gain significantly. 
%
% To prevent this sort of behavior and make Skipping Stone run much faster,
% several computationally costly routines have been translated into
% compiled MATLAB-functions, called MATLAB Executables (MEX). This function 
% will compile all the CPU-intensive and JIT-incompatible routines into a 
% binary executable, suited for your platform. 
%
% Here's a list of all the computationally intesive functions that will be
% compiled by running this function:
%
%   lambert.m             : lambert targeter for ballistic flights
%   lambert_low_exposins.m: lambert targeter for Exponential Sinusoids
%   TOF.m                 : time-of-flight equation for central body flyby's
%   progress_orbit.c      : Kepler state transition matrix
%   eM2E.c                : solve Kepler's equation
%   paretofront.c         : logical Paretofront membership test
%
function speedup

    % Please report bugs and inquiries to:
    %
    % Name       : Rody P.S. Oldenhuis
    % E-mail     : oldenhuis@gmail.com    (personal)
    %              oldenhuis@luxspace.lu  (professional)
    % Affiliation: LuxSpace sarl
    % Licence    : BSD


    % If you find this work useful, please consider a donation:
    % https://www.paypal.me/RodyO/3.5


    %% INITIALIZE

    % codegen was called 'emlmex' before R20103a:
    use_emlmex = verLessThan('MATLAB', '8.1');
        
    % Collector for all warning/error messages; will be reported via
    % warndlg() at the end
    msg = {};

    % get confirmation (compilation is time consuming and uninterruptible)
    while true
        clc, fprintf(1, [...
                     ' This function will attempt to compile all of the following code\n',...
                     ' for the %s platform:\n\n',...
                     '  lambert.m             : lambert targeter for ballistic flights\n',...
                     '  lambert_low_exposins.m: lambert targeter for Exponential Sinusoids\n',...
                     '  TOF.m                 : time-of-flight equation for central body flyby''s\n',...
                     '  progress_orbit.c      : Kepler state transition matrix\n',...
                     '  eM2E.c                : solve Kepler''s equation\n',...
                     '  paretofront.c         : logical Paretofront membership test\n\n',...
                     ' Do you wish to continue? '], computer('arch'));
        yn = input('','s');
        if any(strcmpi(yn, {'y' 'yes' 'yup' 'yep' 'ok' 'yeah' 'jawohl' 'yessirree'})), break;
        elseif any(strcmpi(yn, {'n' 'no' 'nah' 'nope' 'neh' 'nvm' 'nop' 'negative'})), return;
        else continue;
        end
    end

    % make sure the user's MEX-setup is OK
    mex -setup

    % set the current directory equal to where this file is located
    % (normally in the same directory as MAIN.M, and the rest of the
    % Skipping Stone code)
    rootdir = fileparts(mfilename('fullpath'));
    cd(rootdir);

    %% COMPILE ORBITAL MECHANICS ROUTINES

    % Compile Kepler State Transition Matrix routine
    cd([rootdir,filesep,'orbital_mechanics']);
    fprintf(1, 'Compiling PROGRESS_ORBIT.C...\n');
    try
        mex progress_orbit.c
    catch ME,ME; %#ok<VUNUS> (suppresses incorrectly generated warning (#ok<MUCTH>); bug in R2010a)
        msg{end+1} = ['COMPILING PROGRESS_ORBIT.C FAILED: ',...
                      ME.message];
    end

    % compile high-thrust lambert-targeter
    % (Written in Embedded MATLAB, so use ELMEX())
    cd([rootdir,filesep,...
        'orbital_mechanics',filesep,...
        'Lambert targeters',filesep,...
        'high_thrust']);
    
    fprintf(1, 'Compiling LAMBERT.M...\n');
    example_input = {...
        [0.0, 0.0, 0.0], ...% r1vec
        [0.0, 0.0, 0.0], ...% r2vec
        0.0, ...            % tf
        0.0, ...            % m
        0.0};%#ok           % muC
    % compile
    try
        if use_emlmex
            emlmex -eg example_input lambert.m %#ok<DEMLMEX>
        else
            codegen -config:mex -args example_input -o lambert lambert.m
        end
    catch ME,ME; %#ok<VUNUS> 
        msg{end+1} = ['COMPILING LAMBERT.M FAILED: ',...
                      ME.message];
    end

    % compile Low-thrust (ExpoSins) lambert-targeter
    % (Written in Embedded MATLAB, so use ELMEX())
    cd(fullfile(rootdir,...
                'orbital_mechanics',...
                'Lambert targeters',...
                'low thrust'));
    fprintf(1, ...
            'Compiling LAMBERT_LOW_EXPOSINS.M...\n');        
    example_input = {...
        [0.0, 0.0, 0.0],...  % r1vec
        [0.0, 0.0, 0.0],...  % r2vec
        0.0, ...             % tf
        0.0, ...             % k2
        0.0, ...             % N
        0.0, ...             % M0
        0.0, ...             % Isp
        0.0};%#ok            % muC
    
    % compile
    try
        if use_emlmex
            emlmex -eg example_input lambert_low_exposins.m %#ok<DEMLMEX>
        else
            codegen -config:mex -args example_input -o lambert_low_exposins lambert_low_exposins.m
        end
    catch ME,ME; %#ok<VUNUS> 
        msg{end+1} = ['COMPILING LAMBERT_LOW_EXPOSINS.M FAILED: ',...
                      ME.message];
    end

    % compile time-of-flight equation for Central Body Flyby
    % (Written in Embedded MATLAB, so use ELMEX())
    cd(fullfile(rootdir,...
                'orbital_mechanics',...
                'Lambert targeters',...
                'central body flyby'));
    fprintf(1, ...
            'Compiling TOF.M...\n');
    example_input = {...
        0.0,...   % theta
        0.0,...   % pericenter radius
        0.0,...   % second radius
        false,... % compute gradients?
        0.0};%#ok % muC
    
    % compile
    try
        if use_emlmex
            emlmex -eg example_input TOF.m %#ok<DEMLMEX>
        else
            codegen -config:mex -args example_input -o TOF TOF.m
        end
    catch ME,ME; %#ok<VUNUS> 
        msg{end+1} = ['COMPILING TOF.M FAILED: ',...
                      ME.message];
    end

    %% COMPILE CONVERSION ROUTINES

    % Compile conversion routine from mean anomaly to eccentric anomaly
    % (Equal to solving Kepler's equation. Actually, the compiler version
    % is only slightly faster for elliptic orbits; for hyperbolic
    % trajectories, it's surprisingly much slower than MATLAB)
    cd(fullfile(rootdir, 'convert'));
    fprintf(1, 'Compiling eM2E.C...\n');
    % compile
    try
        mex eM2E.c
    catch ME,ME; %#ok<VUNUS> 
        msg{end+1} = ['COMPILING eM2E.C FAILED: ',...
                      ME.message];
    end

    %% COMPILE OPTIMIZATION ROUTINES

    % Compile Pareto front membership tester
    cd(fullfile(rootdir,...
                'optimization',...
                'GODLIKE'));
    fprintf(1,...
            'Compiling PARETOFRONT.C...\n');
    % compile
    try
        mex paretofront.c
    catch ME,ME; %#ok<VUNUS> 
        msg{end+1} = ['COMPILING PARETOFRONT.C FAILED: ',...
                      ME.message];
    end

    %% FINALIZE

    % Display final message
    if isempty(msg)
        fprintf(1,...
                'All functions successfully compiled. Enjoy the speed!\n\n');
    else
        warndlg(msg);
    end

    % goto root dir
    cd(rootdir);

end % speedup
