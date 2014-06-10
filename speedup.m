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
% is not very intelligent; it will for example try to compile nested-loops 
% at each iteration of the outer loop in some cases, or re-compile small 
% loops embedded in functions which also contain code unsuitable for the 
% JIT-compiler at every call of that particular function. This implies 
% significant overhead due to compilation for some functions, and thus also 
% much longer runtimes. 
%
% To prevent this sort of behavior and make Skipping Stone run much faster,
% several computationally costly routines have been translated into 
% compiled MATLAB-functions, called MATLAB Executables (MEX). As with any
% compiled program/library/function, its corresponding compiled version is 
% platform dependent, i.e., if you are running UNIX, routines compiled on
% it will generally NOT run on MacOS, Windows, etc. 
%
% That is where this function comes in. This function will compile all
% the CPU-intensive and JIT-incompatible routines into a format suited for 
% your platform. Currently, all the compilable functions have already been
% compiled to MS Windows 32-bit versions, so if you're running a Windows
% 32-bit operating system, running this function will essentially do
% nothing. 
%
% However, if you're running any other operating system, executing this 
% function once will most likely increase computation speed tremendously.
%
% Here's a list of all the computationally intesive functions that will be
% compiled by running this function:
%
%   lambert_high.m        : lambert targeter for ballistic flights
%   lambert_low_exposins.m: lambert targeter for Exponential Sinusoids
%   TOF.m                 : time-of-flight equation for central body flyby's
%   progress_orbit.c      : Kepler state transition matrix
%   eM2E.c                : solve Kepler's equation
%   paretofront.c         : logical Paretofront membership test 
%
%
%

function speedup
% Last edited 14/Nov/2009

    %% INITIALIZE
    
    % all warning/error messages
    msg = {};
    
    % Check if compilation is neccessary
    if any(strcmpi(computer('arch'), {'Win32'}))
        clc, fprintf(1, 'All routines have already been compiled to %s platform.\n', computer('arch'));
        return
    end
    
    % get confirmation (compilation is lengthy, uninterruptible and annoying)
    while true
        clc, fprintf(1, ...
           [' This function will attempt to compile all of the following code\n',...
            ' for the %s platform:\n\n',...
            '  lambert_high.m        : lambert targeter for ballistic flights\n',...
            '  lambert_low_exposins.m: lambert targeter for Exponential Sinusoids\n',...
            '  TOF.m                 : time-of-flight equation for central body flyby''s\n',...
            '  progress_orbit.c      : Kepler state transition matrix\n',...
            '  eM2E.c                : solve Kepler''s equation\n',...
            '  paretofront.c         : logical Paretofront membership test\n\n',...
            ' Do you wish to continue? '], computer('arch'));
        yn = input('','s');
        if any(strcmpi(yn, {'y';'yes';'yup';'yep';'ok';'yeah'})), break;
        elseif any(strcmpi(yn, {'n';'no';'nope';'neh';'nvm'})), return;
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
    catch ME
        msg{end+1} = ['COMPILING PROGRESS_ORBIT.C FAILED: ' ME.message];
    end
    
    % compile high-thrust lambert-targeter
    % (Written in Embedded MATLAB, so use ELMEX())
    cd([rootdir,filesep,'orbital_mechanics',filesep,'Lambert targeters',filesep,'high thrust']);
    fprintf(1, 'Compiling LAMBERT_HIGH.M...\n');
    example_input = {...
        [0.0, 0.0, 0.0], ...% r1vec
        [0.0, 0.0, 0.0], ...% r2vec
        0.0, ...            % tf
        0.0, ...            % m
        0.0};%#ok           % muC
    % compile
    try
        emlmex -eg example_input lambert_high.m
    catch ME
        msg{end+1} = ['COMPILING LAMBERT_HIGH.M FAILED: ' ME.message];
    end
    
    % compile Low-thrust (ExpoSins) lambert-targeter
    % (Written in Embedded MATLAB, so use ELMEX())
    cd([rootdir,filesep,'orbital_mechanics',filesep,'Lambert targeters',filesep,'low thrust']);
    fprintf(1, 'Compiling LAMBERT_LOW_EXPOSINS.M...\n');
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
        emlmex -eg example_input lambert_low_exposins.m
    catch ME
        msg{end+1} = ['COMPILING LAMBERT_LOW_EXPOSINS.M FAILED: ' ME.message];
    end
    
    % compile time-of-flight equation for Central Body Flyby
    % (Written in Embedded MATLAB, so use ELMEX())
    cd([rootdir,filesep,'orbital_mechanics',filesep,'Lambert targeters',filesep,'central body flyby']);
    fprintf(1, 'Compiling TOF.M...\n');
    example_input = {...
        0.0,...   % theta
        0.0,...   % pericenter radius
        0.0,...   % second radius
        false,... % compute gradients?
        0.0};%#ok % muC
    % compile
    try
        emlmex -eg example_input TOF.m
    catch ME
        msg{end+1} = ['COMPILING TOF.M FAILED: ' ME.message];
    end
    
    %% COMPILE CONVERSION ROUTINES
    
    % Compile conversion routine from mean anomaly to eccentric anomaly
    % (Equal to solving Kepler's equation. Actually, the compiler version 
    % is only slightly faster for elliptic orbits; for hyperbolic 
    % trajectories, it's surprisingly much slower than MATLAB)
    cd([rootdir,filesep,'convert']);
    fprintf(1, 'Compiling eM2E.C...\n');
    % compile
    try
        mex eM2E.c
    catch ME
        msg{end+1} = ['COMPILING eM2E.C FAILED: ' ME.message];
    end
    
    %% COMPILE OPTIMIZATION ROUTINES
    
    % Compile Pareto front membership tester
    cd([rootdir,filesep,'optimization',filesep,'GODLIKE']);
    fprintf(1, 'Compiling PARETOFRONT.C...\n');
    % compile
    try
        mex paretofront.c
    catch ME
        msg{end+1} = ['COMPILING PARETOFRONT.C FAILED: ', ME.message];
    end
    
    %% FINALIZE
    
    % Display final message
    if isempty(msg)
        fprintf(1, 'All functions successfully compiled. Enjoy the speed!');
    else
        warndlg(msg); 
    end
    
    % goto root dir
    cd(rootdir);
    
end % speedup
