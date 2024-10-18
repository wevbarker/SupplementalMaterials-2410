% Get the folder where the script is located
script_folder = fileparts(mfilename('fullpath'));

% List of dataset file names (update this list as needed)
file_names = {'1e8g_GaussianFitted_d_1', '1e8g_GaussianFitted', '1e8g_NewFitted_n_4', '1e8g_NewFitted'}; % Add your file names here

% Frequency region to plot GW signature
lowF = 1e-2;
highF = 3e6;
conv = 1.546e-15;
k1 = logspace(log10(lowF / conv), log10(highF / conv), 250);

% Constants
cg = 0.4;
OmegaR0 = 2.47e-5;

% Loop over all provided file names
for i = 1:length(file_names)
    % Build the filename for the power spectrum and model data
    variant_name = file_names{i};
    power_spectrum_file = fullfile(script_folder, [variant_name, '.txt']);
    model_data_file = fullfile(script_folder, [variant_name, '_GW.csv']);
    
    % Load the power spectrum data
    data = load(power_spectrum_file);
    k_list = data(:, 1);
    spectrum = data(:, 2);
    
    % Interpolate the power spectrum
    logPInterpolation = @(x) interp1(log10(k_list), log10(spectrum), log10(x), 'linear', 'extrap');
    P = @(k) 10.^logPInterpolation(k);

    % Functions Ic and Is
    Ic = @(d, s) -36 * pi * (((s.^2 + d.^2 - 2).^2 ./ (s.^2 - d.^2).^3) .* (s >= 1));
    Is = @(d, s) -36 * ((s.^2 + d.^2 - 2) ./ (s.^2 - d.^2).^2) .* ...
        ((s.^2 + d.^2 - 2) ./ (s.^2 - d.^2) .* log((1 - d.^2) ./ abs(s.^2 - 1)) + 2);

    % Integrand for GW signature
    Integrand2 = @(k, d, s) ((d.^2 - 1/3) .* (s.^2 - 1/3) ./ (s.^2 - d.^2)).^2 .* ...
        P(k * sqrt(3) * (s + d) / 2) .* P(k * sqrt(3) * (s - d) / 2) .* ...
        (Is(d, s).^2 + Ic(d, s).^2);

    % GW signature calculation
    h2OmegaGW2 = @(k) integral2(@(d, s) Integrand2(k, d, s), 0, 1/sqrt(3), 1/sqrt(3), Inf, ...
        'AbsTol', 1e-23, 'Method', 'iterated') * cg * OmegaR0 / 36;

    % Calculate h2OmegaGW2 for each k1
    h2OmegaGW2Sol = arrayfun(@(k) h2OmegaGW2(k), k1);

    % Combine the frequency and corresponding GW signature into a single matrix
    output_data = [conv * k1', h2OmegaGW2Sol'];

    % Save the data to a CSV file in the same folder as the script
    output_csv_file = fullfile(script_folder, [variant_name, '_GW_Matlab.csv']);
    writematrix(output_data, output_csv_file);
end
