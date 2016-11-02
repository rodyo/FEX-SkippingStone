function [model, constants] = minor_planets_model(model, constants, environment, reload)
%MINOR_PLANET_MODEL   Definition file that loads the compete model for
%                     the Solar system's minor planets

% Author:
% Name       : Rody P.S. Oldenhuis
% E-mail     : oldenhuis@gmail.com

% If you find this work useful, please consider a small donation:
% https://www.paypal.me/RodyO/3.5


    % keep track of whether the model has been loaded already
    persistent loaded
    if isempty(loaded)
        loaded = false; end

    % set default re-load status
    if nargin == 3
        reload = false; end

    progress_prefix = 'Loading MP model; ';

    % only load if not loaded already
    if ~loaded || reload

        % initialize model
        %{
        zero  = zeros(model.MPs.number_of_MPs, 1);
        zero2 = cell (model.MPs.number_of_MPs, 1);

        model.MPs.types          = zero;  model.MPs.es     = zero;
        model.MPs.M0s            = zero;  model.MPs.omegas = zero;
        model.MPs.Omegas         = zero;  model.MPs.is     = zero;
        model.MPs.as             = zero;  model.MPs.epochs = zero;
        model.MPs.numbers        = zero2;
        model.MPs.first_observed = zero;  model.MPs.last_observed          = zero;
        model.MPs.oppositions    = zero;  model.MPs.number_of_observations = zero;
        model.MPs.epoch_MJD      = zero;  model.MPs.epoch_days_past_J2000  = zero;
        %}


        % EXTRACT DATA FROM DATABASE
        % ======================================================================

        % initialize progress bar
        progress_bar(0, [progress_prefix 'reading asteroid data...']);
        pause(0.25); % allow some time to update the progress bar

        % initialize loop
        fid         = fopen(environment.pathing.MP_filename);
        OC          = onCleanup(@() any(fopen('all')==fid) && fclose(fid) );
        headerlines = 41;
        numcolumns  = 202;

        % Estaimte total amount of Minor Planets
        %total = countlines(filename) - headerlines;
        % Since the format is so predictable, the quickest way to
        % estimate the number of MPs in the file is:
        num_MPs = dir(environment.pathing.MP_filename);
        num_MPs = ceil( num_MPs.bytes/numcolumns );

        % Initialize string data container
        astdata(num_MPs, numcolumns) = char(0);

        % Read file in batches of 25.000 lines to prevent out-of-memory
        % errors, allowing us to display progress as well
        counter   = 0;
        batch_MPs = 2.5e4;
        while ~feof(fid)

            % skip the header in the first read
            new_MPs = textscan(fid,...
                               '%s',...
                               batch_MPs, ...
                               'headerlines', headerlines * (counter == 0),...
                               'delimiter'  , '\n',...
                               'bufsize'    , 65535,...
                               'Whitespace' , '');

            % Strip empty/newlines
            new_MPs = strtrim(new_MPs{1});
            new_MPs = new_MPs(~cellfun('isempty', new_MPs));

            astdata( (1:numel(new_MPs)) + counter, :) = char(new_MPs);

            counter = counter + numel(new_MPs);
            progress_bar(0.5 * counter/num_MPs, [progress_prefix,...
                         num2str(counter), ' read...']);
        end
        fclose(fid);

        % Save accurate count
        model.MPs.number_of_MPs = counter;

        % Chop off insignificant trailing zeros
        astdata(counter+1:end,:) = [];

        % Reading done; commence the conversion
        progress_bar(0.5, [progress_prefix,...
                     'read complete; ' num2str(counter) ' minor planets in file.']);

        % extract all relevant info
        progress_bar(0.25, [progress_prefix 'converting asteroid data...']);

        % NOTE that the "magic numbers" appearing in the extraction, are
        % documented in the header of the minor planets datafile; They
        % indicate the column numbers the specific data can be located.
        numbers = astdata(:,   1:  7);  % asteroid numbers (contains numbers+alpha-numeric)
        epochs  = astdata(:,  21: 25);  % epochs
        M0s     = astdata(:,  27: 35);  % mean anomalies
        omegas  = astdata(:,  38: 46);  % arguments of Perihelion
        Omegas  = astdata(:,  49: 57);  % longitudes of Ascending Nodes
        is      = astdata(:,  60: 68);  % inclinations
        es      = astdata(:,  71: 79);  % eccentricities
        as      = astdata(:,  93:103);  % semi-major axes
        types   = astdata(:, 162:165);  % Asteroid type

        % The names can be inerted as-is
        model.MPs.MP_names = astdata(:, 167:193);

        oppositions            = astdata(:, 124:126);  % number of oppositions
        first_observed         = astdata(:, 128:131);  % year of first observation
        last_observed          = astdata(:, 133:136);  % year of last observation
        number_of_observations = astdata(:, 118:122);  % number of observations

        clear astdata

        % oOoOoOoOoOoOoOoOoOoOoOoOoOoOoOoOoOoOoOoOoOoOoOoOoOoOoOoOoOoOoOoOoOoOoOoO
        % FROM THE MPCORB-DOCUMENTATION:
        % oOoOoOoOoOoOoOoOoOoOoOoOoOoOoOoOoOoOoOoOoOoOoOoOoOoOoOoOoOoOoOoOoOoOoOoO
        %
        % Dates of the form YYYYMMDD may be packed into five characters to conserve
        % space.
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
        century = 100*(epochs(:,1) - 55);
        % Complete the year
        year = century + quick_str2num(epochs(:,2:3));

        % find the correct month and day
        % asically, '1'-'9' map to 1-9, and 'A'-'V' map to 10-31:
        md = epochs(:,4:5) - '0';
        md(md>9) = md(md>9) - 'A' + '0' + 10;

        month = md(:,1);
        day   = md(:,2);

        clear epochs

        % convert the time
        MJD = date2MJD(year, month, day);
        days_past_J2000 = date2days(year, month, day);

        % insert in model
        model.MPs.epoch_MJD = MJD;
        model.MPs.epoch_days_past_J2000 = days_past_J2000;

        % convert everything else
        d2r = pi/180;
        model.MPs.numbers = numbers;
        model.MPs.M0s     = quick_str2num(M0s)    * d2r; % convert to [rad]
        model.MPs.omegas  = quick_str2num(omegas) * d2r;
        model.MPs.Omegas  = quick_str2num(Omegas) * d2r;
        model.MPs.is      = quick_str2num(is)     * d2r;
        model.MPs.es      = quick_str2num(es);
        model.MPs.as      = quick_str2num(as)     * constants.AU; % convert to [km]
        model.MPs.types   = hex2dec(types);                 % convert to decimals
        model.MPs.oppositions            = quick_str2num(oppositions);
        model.MPs.number_of_observations = quick_str2num(number_of_observations);

        % sometimes, especially for the newer additions, the arc is
        % smaller than one year. In those cases, the arclength is given
        % as a number of days. This translates into:
        model.MPs.first_observed = quick_str2num(first_observed);
        model.MPs.last_observed  = quick_str2num(last_observed);
        % ('days' converts to 57697)
        model.MPs.last_observed(model.MPs.last_observed > 52e3) = 0;

        % All string conversions done; these big honking piles of strings
        % are no longer needed
        clear numbers M0s omegas Omegas is es as oppositions
        clear number_of_observations last_observed first_observed


        % CALCULATIONS
        % =================================================================

        progress_bar(0.75, [progress_prefix 'calculating initial values...']);

        % true anomalies
        model.MPs.theta0s = eM2theta(model.MPs.es, model.MPs.M0s);
        % eccentric anomalies
        model.MPs.E0s = eM2E(model.MPs.es, model.MPs.M0s);
        % x0s
        model.MPs.x0s = kep2cart([model.MPs.as,...
                                  model.MPs.es,...
                                  model.MPs.is,...
                                  model.MPs.Omegas, ...
                                  model.MPs.omegas,...
                                  model.MPs.M0s], constants.mu_central);

        progress_bar(0.85, [progress_prefix 'calculating initial values...']);

        % mean motions
        model.MPs.ns = sqrt(constants.mu_central ./ model.MPs.as.^3);
        % orbital periods
        model.MPs.Ts = 2*pi ./ model.MPs.ns;
        % angular momentum vectors
        model.MPs.Hs = cross( model.MPs.x0s(:, 1:3), model.MPs.x0s(:, 4:6));
        % magnitudes angular momentum vectors
        model.MPs.Hnorms = sqrt(sum(model.MPs.Hs.^2, 2));
        % apses
        model.MPs.peris = model.MPs.as.*(1 - model.MPs.es);
        model.MPs.apos  = model.MPs.as.*(1 + model.MPs.es);
        % parametrization
        progress_bar(0.95, [progress_prefix 'calculating initial values...']);
        [ignored, model.MPs.bs, model.MPs.Cs, model.MPs.Us, model.MPs.Vs] = ...
        kep2para([model.MPs.as, model.MPs.es, model.MPs.is, model.MPs.Omegas, ...
                  model.MPs.omegas, model.MPs.theta0s], constants.mu_central, 'theta');%#ok

        % All done
        progress_bar(1, [progress_prefix ,...
                     'Data for ', num2str(model.MPs.number_of_MPs),' loaded into memory.']);
        pause(1)
        progress_bar('');
        loaded = true;
    end

end % minor planets model

% SUPERFAST conversion from string to double
% (It can be this fast because it isn't very thorough)
function num = quick_str2num(str)

    % Convert to numeric form
    nmstr = str - '0';

    % nullify all spaces
    nmstr(nmstr == 32-'0') = 0;

    % Find position of the decimal points
    pt = find(str(1,:) == '.', 1, 'first');

    % No point; numbers are ints
    if isempty(pt)
        num = nmstr * 10.^(size(str,2)-1 : -1 : 0)';

    % Floats
    else
        % Create powers of ten, padded with a zero at the decimal point
        pows = [10.^(pt-2:-1:0) 0 10.^(-1:-1:-(size(str,2)-pt))].';

        % Check if all periods are in the same place
        % If so, the conversion is fastest
        pt_OK = str(:,pt) == '.';
        if all(pt_OK)
            num = nmstr * pows;

        % But abberant strings need recursion
        else
            num = zeros(size(nmstr,1),1);
            num( pt_OK) = nmstr(pt_OK,:) * pows;
            num(~pt_OK) = quick_str2num(str(~pt_OK,:));
        end
    end

end
