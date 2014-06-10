function [model, constants] = minor_planets_model(model, constants, environment, reload)
%MINOR_PLANET_MODEL   Definition file that loads the compete model for
%                     the Solar system's minor planets
    
% Authors
% .·`·.·`·.·`·.·`·.·`·.·`·.·`·.·`·.·`·.·`·.·`·.·`·.·`·.·`·.·`·.
% Name       : Rody P.S. Oldenhuis
% E-mail     : oldenhuis@dds.nl / oldenhuis@gmail.com
% Affiliation: Delft University of Technology
%              Partly performed at ISAS/JAXA
% .·`·.·`·.·`·.·`·.·`·.·`·.·`·.·`·.·`·.·`·.·`·.·`·.·`·.·`·.·`·.
    
    % Last Edited: 23/Sep/2009
    
    % keep track of whether the model has been loaded already
    persistent loaded
    if isempty(loaded), loaded = false; end
   
    % set default re-load status
    if nargin == 3, reload = false; end
        
    % helper strings
    month_help = {'1';'2';'3';'4';'5';'6';'7';'8';'9';'A';'B';'C'};
    month_nums = 1:12;
    day_help   = {'1';'2';'3';'4';'5';'6';'7';'8';'9';'A';'B';'C';
        'D';'E';'F';'G';'H';'I';'J';'K';'L';'M';'N';'O';'P';'Q';'R';'S';'T';'U';'V'};            
    day_nums = 1:31;
    
    % only load if not loaded already
    if ~loaded || reload
                
        % initialize model
        zero  = zeros(model.MPs.number_of_MPs, 1); 
        zero2 = cell(model.MPs.number_of_MPs, 1); 
        model.MPs.types          = zero; model.MPs.es     = zero;           
        model.MPs.M0s            = zero; model.MPs.omegas = zero;            
        model.MPs.Omegas         = zero; model.MPs.is     = zero;          
        model.MPs.as             = zero; model.MPs.epochs = zero;    
        model.MPs.numbers        = zero2;         
        model.MPs.first_observed = zero; model.MPs.last_observed          = zero;
        model.MPs.oppositions    = zero; model.MPs.number_of_observations = zero;
        model.MPs.epoch_MJD      = zero; model.MPs.epoch_days_past_J2000  = zero;
        
        % XxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXx
        % EXTRACT DATA FROM DATABASE 
        % XxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXx
        
        % initialize progress bar
        progress_bar(0, 'Reading & converting asteroid data...');
        pause(0.25); % allow some time to update the progress bar
                        
        % initialize loop
        fid         = fopen(environment.pathing.MP_filename); 
        index       = 0; 
        total       = model.MPs.number_of_MPs;  
        headerlines = 39;
        
        % first skip the header
        for i = 1:headerlines, fgetl(fid); end
        
        % then loop through the file, one line at a time
        % NOTE: this is the only way ALL data gets extracted. The canned
        % routine TEXTSCAN() will produce errors on empty lines, and some 
        % lines contain data in columns where there is none in other lines,
        % also causing TEXTSCAN() to error out. FGETL() is the only thing
        % that does the trick.
        while ~feof(fid)
            % get one line
            astdata = fgetl(fid);               
            % some lines are empty, so skip these
            if isempty(astdata), continue, end            
            % increase index counter 
            index = index + 1; 
            
            % extract all relevant info
            
            % NOTE that the "magic numbers" appearing in the extraction, are
            % documented in the header of the minor planets datafile; They
            % indicate the column numbers the specific data can be located.
            numbers = astdata(1:7);             % asteroid numbers (contains numbers+alpha-numeric)
            epochs  = astdata(21:25);           % epochs
            M0s     = astdata(27:35);           % mean anomalies            
            omegas  = astdata(38:46);           % arguments of Perihelion            
            Omegas  = astdata(49:57);           % longitudes of Ascending Nodes            
            is      = astdata(60:68);           % inclinations            
            es      = astdata(71:79);           % eccentricities            
            as      = astdata(93:103);          % semi-major axes            
            number_of_observations = astdata(118:122);   % number of observations
            types   = astdata(162:165);         % Asteroid type                    
            oppositions    = astdata(124:126);  % number of oppositions                   
            first_observed = astdata(128:131);  % year of first observation  
            last_observed  = astdata(133:136);  % year of last observation
            
            % oOoOoOoOoOoOoOoOoOoOoOoOoOoOoOoOoOoOoOoOoOoOoOoOoOoOoOoOoOoOoOoOoOoOoOoO
            % FROM THE MPCORB-DOCUMENTATION:
            % oOoOoOoOoOoOoOoOoOoOoOoOoOoOoOoOoOoOoOoOoOoOoOoOoOoOoOoOoOoOoOoOoOoOoOoO
            %
            % Dates of the form YYYYMMDD may be packed into five characters to conserve 
            % space.
            
            % NOTE: STORING MJD'S FOR THE APPLICABLE RANGES WOULD IMPLY A WHOPPING 0.0053% 
            % INCREASE IN FILE SIZE. GREAT FEATURE, THANKS FOR THAT!
            
            %
            % The first two digits of the year are packed into a single character in column 1
            % (I = 18, J = 19, K = 20). Columns 2-3 contain the last two digits of the year.
            % Column 4 contains the month and column 5 contains the day, coded as detailed below:
            %
            %    Month     Day      Character         Day      Character
            %                      in Col 4 or 5              in Col 4 or 5
            %    Jan.       1           1             17           H
            %    Feb.       2           2             18           I
            %    Mar.       3           3             19           J
            %    Apr.       4           4             20           K
            %    May        5           5             21           L
            %    June       6           6             22           M
            %    July       7           7             23           N
            %    Aug.       8           8             24           O
            %    Sept.      9           9             25           P
            %    Oct.      10           A             26           Q
            %    Nov.      11           B             27           R
            %    Dec.      12           C             28           S
            %              13           D             29           T
            %              14           E             30           U
            %              15           F             31           V
            %              16           G
            %
            % Examples:
            %
            %    1996 Jan. 1    = J9611
            %    1996 Jan. 10   = J961A
            %    1996 Sept.30   = J969U
            %    1996 Oct. 1    = J96A1
            %    2001 Oct. 22   = K01AM
            %
            % This system can be extended to dates with non-integral days. The decimal fraction
            % of the day is simply appended to the five characters defined above.
            %
            % Examples:
            %
            %    1998 Jan. 18.73     = J981I73
            %    2001 Oct. 22.138303 = K01AM138303
            
            % trick to find the century
            century = 100*(epochs(1) - 55);
            % find the year
            year = sscanf(epochs(2:3),'%f',1);
            % complete the year
            year = century + year;
            % find the correct month and day
            month = month_nums(strcmpi(epochs(4), month_help));
            day = day_nums(strcmpi(epochs(5), day_help));
            % convert the time
            MJD = date2MJD(year, month, day, 0, 0, 0);
            days_past_J2000 = date2days(year, month, day, 0, 0, 0);
            % insert in model
            model.MPs.epoch_MJD(index) = MJD;
            model.MPs.epoch_days_past_J2000(index) = days_past_J2000;
                
            % convert everything else
            model.MPs.numbers{index}= numbers;
            model.MPs.M0s(index)    = sscanf(M0s,'%f',1)    * pi/180;       % convert to [rad]            
            model.MPs.omegas(index) = sscanf(omegas,'%f',1) * pi/180;       % convert to [rad]            
            model.MPs.Omegas(index) = sscanf(Omegas,'%f',1) * pi/180;       % convert to [rad]            
            model.MPs.is(index)     = sscanf(is,'%f',1)     * pi/180;       % convert to [rad]            
            model.MPs.es(index)     = sscanf(es,'%f',1);            
            model.MPs.as(index)     = sscanf(as,'%f',1)     * constants.AU; % convert to [km]            
            model.MPs.types(index)  = hex2dec(types);                       % convert to decimals
            model.MPs.oppositions(index)    = sscanf(oppositions,'%f',1);
            model.MPs.number_of_observations(index) = sscanf(number_of_observations,'%f',1);
            
            % sometimes, especially for the newer additions, the arc is
            % smaller than one year. In those cases, the arclength is given 
            % as a number of days. This translates into:
            model.MPs.first_observed(index) = sscanf(first_observed,'%f',1); 
            lo = sscanf(last_observed,'%f',1);
            if isempty(lo)
                model.MPs.last_observed(index) = 0;
            else
                model.MPs.last_observed(index) = lo;
            end           
            
            % progress bar
            if(mod(index, 1e4) == 0)
                % max 99%
                progress_bar(index/total-0.01, 'Reading & converting asteroid data...');
            end                        
        end        
                
        % XxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXx
        % CALCULATIONS
        % XxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXxXx
                
        % initial thetas
        progress_bar(0.99, 'Calculating initial values...');        
        % true anomalies
        model.MPs.theta0s = eM2theta(model.MPs.es, model.MPs.M0s);
        % eccentric anomalies        
        model.MPs.E0s = eM2E(model.MPs.es, model.MPs.M0s);
        % calculate x0s        
        model.MPs.x0s = kep2cart([model.MPs.as, model.MPs.es, model.MPs.is, model.MPs.Omegas, ...
                       model.MPs.omegas, model.MPs.M0s], constants.mu_central);        
        % calculate mean motions        
        model.MPs.ns = sqrt(constants.mu_central ./ model.MPs.as.^3);        
        % calculate orbital periods        
        model.MPs.Ts = 2*pi ./ model.MPs.ns;        
        % calculate angular momentum vectors        
        model.MPs.Hs = cross( model.MPs.x0s(:, 1:3), model.MPs.x0s(:, 4:6));        
        % calculate magnitudes angular momentum vectors        
        model.MPs.Hnorms = sqrt(sum(model.MPs.Hs.^2, 2));        
        % calculate apses        
        model.MPs.peris = model.MPs.as.*(1 - model.MPs.es);
        model.MPs.apos  = model.MPs.as.*(1 + model.MPs.es);        
        % parametrization        
        [ignored, model.MPs.bs, model.MPs.Cs, model.MPs.Us, model.MPs.Vs] = ...
        kep2para([model.MPs.as, model.MPs.es, model.MPs.is, model.MPs.Omegas, ...
                  model.MPs.omegas, model.MPs.theta0s], constants.mu_central, 'theta');%#ok
                 
        % close file
        fclose(fid);
        
        % all done
        progress_bar(1, ['All data on ',num2str(model.MPs.number_of_MPs),' loaded into memory.']);
        pause(1), progress_bar('');
        
        % set loaded state
        loaded = true;
        
    end
    
end % minor planets model
