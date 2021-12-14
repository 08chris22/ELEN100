%% ELEN 100L (Electric Circuits II): Project 1, Christian Garcia and Logan Barnes
%
%
% *To Submit:*
%
% # Scanned copy of ALL hand calculations.
% # A MATLAB script and publish the solution using MATLAB's *publish*
% feature.%
% # Turn in all MATLAB files.
% # Turn in all LTSpice files.
% # Turn in Screenshot of LTSpice measurement.
%

%% Initialize MATLAB Environment
%

clear; clc; clf; cla; close all;
format long; format compact;

%% Setup global variables
%

% These Ideal Design element values are fixed in the circuit.
VG  =  1;                   % Generator voltage
R1_ideal =      1250          ;    % Ohms
R2_ideal =      1333          ;    % Ohms
C1_ideal =      0.1e-6          ;    % Farads
C2_ideal =      0.1e-6          ;    % Farads

% These Actual Design element values are fixed in the circuit.
R1_actual =    1300       ;   % Ohms
R2_actual =    1300       ;   % Ohms
C1_actual =    0.1e-6       ;   % Farads
C2_actual =    0.1e-6       ;   % Farads

% Setup values for the poles.
w0 =     3000   ;        % Radians/Second
w1 =     20000   ;        % Radians/Second
f0 =     w0/(2*pi)   ;        % Hertz
f1 =     w1/(2*pi)   ;        % Hertz

% Build an array for the angular frequency and convert it to Hertz.
dw = 10;                   % Step size for analysis
w = [1:dw:w0-dw, ...
     w0, ...
     w0+dw:dw:w1-dw, ...
     w1, ...
     w1+dw:dw:1.0e6];      % Radians/Second (ensure poles are included)

 f =      w/(2*pi)     ;            % Hertz

% These values are used for plotting purposes.
fignum = 1;
plot_left   = 1;       plot_right = 2e5;    % x-axis range (Hertz)
plot_bottom = -90;     plot_top   = 5;      % y-axis range (dB)

%% Problem 3

% Display the component values for the Ideal and Actual designs.
%

disp(' ');
disp('The Ideal Design component values are:');
fprintf('    R1 = %+11.4f Ohms.\n',  R1_ideal  );
fprintf('    R2 = %+11.4f Ohms.\n',  R2_ideal  );
fprintf('    C1 = %+11.4e Farads.\n', C1_ideal  );
fprintf('    C2 = %+11.4e Farads.\n', C2_ideal  );

disp(' ');
disp('The Actual Design component values are:');
fprintf('    R1 = %+11.4f Ohms.\n',  R1_actual  );
fprintf('    R2 = %+11.4f Ohms.\n',  R2_actual  );
fprintf('    C1 = %+11.4e Farads.\n', C1_actual );
fprintf('    C2 = %+11.4e Farads.\n', C2_actual );

%%
% Compute the percent differences between the Ideal and Actual design
% component values.
%

diff_R1_ideal_actual = ( R1_actual - R1_ideal )/abs(R1_ideal)*100;
diff_R2_ideal_actual = ( R2_actual - R2_ideal )/abs(R2_ideal)*100;
diff_C1_ideal_actual = ( C1_actual - C1_ideal )/abs(C1_ideal)*100;
diff_C2_ideal_actual = ( C2_actual - C2_ideal )/abs(C2_ideal)*100;

disp(' ');
disp('The percent difference between Ideal and Actual design');
disp('component values:');
fprintf('    %% diff R1 = %+8.4f (%%).\n', diff_R1_ideal_actual );
fprintf('    %% diff R2 = %+8.4f (%%).\n', diff_R2_ideal_actual );
fprintf('    %% diff C1 = %+8.4f (%%).\n', diff_C1_ideal_actual );
fprintf('    %% diff C2 = %+8.4f (%%).\n', diff_C2_ideal_actual );

%%
% Display the poles for the target circuit design transfer function.
%

disp(' ');
disp('The poles for the circuit are:');
fprintf('    w0 = %+11.4f Radians/Second.\n',w0);
fprintf('    w1 = %+11.4f Radians/Second.\n',w1);
fprintf('    f0 = %+11.4f Hertz.\n',f0);
fprintf('    f1 = %+11.4f Hertz.\n',f1);

%%
% Setup the matrices used to generate the Bode plots for the Ideal and
% Actual designs.
%

G1_ideal = [ ...
      (1)            (0)                         (0); ...
      (-1/R1_ideal)  (1/R1_ideal + 1/R2_ideal)   (-1/R2_ideal); ...
      (0)            (-1/R2_ideal)               (1/R2_ideal)];

G3_ideal = [ 0,0,0;0,0,0;0,0,0 ];

G2_ideal = [ 0,0,0;0,C1_ideal,0;0,0,C2_ideal ];

G1_actual = [ ...
      (1)            (0)                         (0); ...
      (-1/R1_actual)  (1/R1_actual + 1/R2_actual)   (-1/R2_actual); ...
      (0)            (-1/R2_actual)               (1/R2_actual) ];

G3_actual = [ 0,0,0;0,0,0;0,0,0 ];

G2_actual = [ 0,0,0;0,C1_actual,0;0,0,C2_actual ];

B = [ 1;0;0 ];

%%
% Locate the poles in the frequency vector for plotting purposes.
%

% Find the pole values.
pole_1 = 0;
for iter = 1:length(f)             % Locate the first pole
    if (f(iter) == f0)
        pole_1 = iter;
        break;
    end
end

pole_2 = 0;
for iter = pole_1+1:length(f)     % Locate the second pole
    if (f(iter) == f1)
        pole_2 = iter;
        break;
    end
end

%%
% Calculate the frequency response for the Ideal and Actual designs.
%

Hw_ideal   = proj1E100_freqresp( G1_ideal,G2_ideal,G3_ideal,B,w,VG );
Hw_actual  = proj1E100_freqresp( G1_actual,G2_actual,G3_actual,B,w,VG );

% Capture the values at the poles.
Hw_ideal_f0  = Hw_ideal(pole_1);
Hw_ideal_f1  = Hw_ideal(pole_2);

Hw_actual_f0 = Hw_actual(pole_1);
Hw_actual_f1 = Hw_actual(pole_2);

%%
% Generate the plot for $\displaystyle{H_{ideal}(f)}$ and indicate where
% the two poles occur.

fignum = fignum+1; figObj = figure(fignum);  % Establish a figure number
set(fignum, 'Name',['H(f) Ideal Design']);   % Name the figure

Hw_ideal_Plot = semilogx( f, Hw_ideal ,'-r');          % Generate plot
grid on;                                     % Turn grid on
xlabel('Frequency (Hz)');                    % Label the x-axis
ylabel('Amplitude (dB)');                    % Label the y-axis
axis([plot_left, plot_right, ...
      plot_bottom, plot_top]);               % Bound plot
title(['Figure ',num2str(fignum,'%-2.u'),...
       ': H_i_d_e_a_l(f)']);
legend('H_i_d_e_a_l(f)', 'Location', 'NorthEast');

% Add cursors to the plot.
makedatatip(Hw_ideal_Plot, [pole_1; pole_2]);

%%
% Generate the plot for $\displaystyle{H_{actual}(f)}$ and indicate where
% the two poles occur.

fignum = fignum+1; figObj = figure(fignum);  % Establish a figure number
set(fignum, 'Name',['H(f) Actual Design']);  % Name the figure

Hw_actual_Plot = semilogx( f , Hw_actual ,'-b');         % Generate plot
grid on;                                     % Turn grid on
xlabel('Frequency (Hz)');                    % Label the x-axis
ylabel('Amplitude (dB)');                    % Label the y-axis
axis([plot_left, plot_right, ...
      plot_bottom, plot_top]);               % Bound plot
title(['Figure ',num2str(fignum,'%-2.u'),...
       ': H_a_c_t_u_a_l(f)']);
legend('H_a_c_t_u_a_l(f)', 'Location', 'NorthEast');

% Add cursors to the plot.
makedatatip(Hw_actual_Plot, [pole_1; pole_2]);

%%
% Generate the plot for comparing $\displaystyle{H_{ideal}(f)}$ and
% $\displaystyle{H_{actual}(f)}$.

fignum = fignum+1; figObj = figure(fignum);  % Establish a figure number
set(fignum, 'Name', ...
    ['H(f) Ideal and Actual Design']);       % Name the figure

Hw_ideal_actual_Plot = ...
    semilogx( f , Hw_ideal ,'-r', ...
              f , Hw_actual ,'-b');                      % Generate plot
grid on;                                     % Turn grid on
xlabel('Frequency (Hz)');                    % Label the x-axis
ylabel('Amplitude (dB)');                    % Label the y-axis
axis([plot_left, plot_right, ...
      plot_bottom, plot_top]);               % Bound plot
title(['Figure ',num2str(fignum,'%-2.u'),...
       ': H_i_d_e_a_l(f) and H_a_c_t_u_a_l(f)']);
legend('H_i_d_e_a_l(f)', 'H_a_c_t_u_a_l(f)', 'Location', 'NorthEast');

%%
% Calculate the percent difference between $\displaystyle{H_{ideal}(f)}$
% and $\displaystyle{H_{actual}(f)}$ at the two poles.
%

diff_0_ideal_actual = ( Hw_actual_f0 - Hw_ideal_f0 )/abs(Hw_ideal_f0)*100;
diff_1_ideal_actual = ( Hw_actual_f1 - Hw_ideal_f1 )/abs(Hw_ideal_f1)*100;

disp(' ');
disp('The difference between Ideal and Actual designs at the poles:');
fprintf('    Ideal  Design H(%+10.4f) = %+8.4f (dB).\n', f0, Hw_ideal_f0);
fprintf('    Actual Design H(%+10.4f) = %+8.4f (dB).\n', f0, Hw_actual_f0);
fprintf('        %% diff = %+8.4f (%%).\n', diff_0_ideal_actual);
fprintf('    Ideal  Design H(%+10.4f) = %+8.4f (dB).\n', f1, Hw_ideal_f1 );
fprintf('    Actual Design H(%+10.4f) = %+8.4f (dB).\n', f1, Hw_actual_f1 );
fprintf('        %% diff = %+8.4f (%%).\n', diff_1_ideal_actual );

%% Problem 4

% Calculate the percent difference between $\displaystyle{H_{actual}(f)}$
% and $\displaystyle{H_{LTSpice}(f)}$ actual designs at the two poles.

Hw_ltspice_f0 =     -3.167      ;      % dB
Hw_ltspice_f1 =     -19.73      ;      % dB
f0_ltspice =      3000/(2*pi)    ;          % Hertz
f1_ltspice =      20000/(2*pi)    ;          % Hertz

diff_0_actual_ltspice =  ( Hw_ideal_f0 - Hw_actual_f0 )/abs(Hw_actual_f0)*100;
diff_1_actual_ltspice =  ( Hw_ideal_f1 - Hw_actual_f1 )/abs(Hw_actual_f1)*100;

disp(' '); 
disp('The percent difference between MATLAB and LTSpice Actual');
disp('designs at the poles:');
fprintf('    Actual  MATLAB H(%+10.4f) = %+8.4f (dB).\n', ....
         f0, Hw_actual_f0 );
fprintf('    Actual LTSpice H(%+10.4f) = %+8.4f (dB).\n', ...
         f0_ltspice, Hw_ltspice_f0);
fprintf('        %% diff = %+8.4f (%%).\n', diff_0_actual_ltspice);
fprintf('    Actual  MATLAB H(%+10.4f) = %+8.4f (dB).\n', ...
         f1 , Hw_actual_f1 );
fprintf('    Actual LTSpice H(%+10.4f) = %+8.4f (dB).\n', ...
         f1_ltspice , Hw_ltspice_f1 );
fprintf('        %% diff = %+8.4f (%%).\n', diff_1_actual_ltspice );

%% Problem 5
%
% Vary the Actual design component values and calculate the frequency
% response for each variation.
%

% Declare the number of component value iterations.
value_sets =  5 ;

% Build the actual component vector
Q_actual = [R1_actual, R2_actual, C1_actual, C2_actual];

% Generate the frequency response values for the specified number of
% iterations.
[Hw_actual_varied, Q_actual_varied] = ...
             proj1E100_freqresp_varied( Q_actual,B,w,VG,value_sets );

% Capture the values at the poles.
Hw_actual_varied_f0 = zeros(1,value_sets);
Hw_actual_varied_f1 = zeros(1,value_sets);
for iter = 1:value_sets
    Hw_actual_varied_f0(iter)  = Hw_actual_varied(iter, pole_1);
    Hw_actual_varied_f1(iter)  = Hw_actual_varied(iter, pole_2);
end

%%
% Generate the plot for variations in the Actual design component values
% and display all $\displaystyle{H_{varied}(f)}$ curves on a single plot.

fignum = fignum+1; figObj = figure(fignum);  % Establish a figure number
set(fignum, 'Name', ...
    ['H(f) Actual Design Varied']);          % Name the figure

Hw_actual_varied_Plot = ...
    semilogx( f , Hw_actual_varied );            % Generate plot
grid on;                                     % Turn grid on
xlabel('Frequency (Hz)');                    % Label the x-axis
ylabel('Amplitude (dB)');                    % Label the y-axis
axis([plot_left, plot_right, ...
      plot_bottom, plot_top]);               % Bound plot
title(['Figure ',num2str(fignum,'%-2.u'),...
       ': Varied H_a_c_t_u_a_l(f)']);
legend('H_1(f)', 'H_2(f)', 'H_3(f)', 'H_4(f)', 'H_5(f)', ...
       'Location', 'NorthEast');
   
%%
% Calculate the percent difference between $\displaystyle{H_{varied}(f)}$
% and $\displaystyle{H_{actual}(f)}$ at the two poles of each variation.
%

diff_0_actual_varied = ...
    (Hw_actual_varied_f0-Hw_actual_f0)/abs(Hw_actual_f0)*100;
diff_1_actual_varied = ...
    ( Hw_actual_varied_f1-Hw_actual_f1 )/abs( Hw_actual_f1 )*100;

disp(' ');
disp('The difference between Varied and Actual designs at the poles:');
for iter = 1:value_sets
    diff_R1_actual_varied = ...
        (Q_actual_varied(iter,1)-R1_actual)/abs(R1_actual)*100;
    diff_R2_actual_varied = ...
        (Q_actual_varied(iter,2)- R2_actual )/abs( R2_actual )*100;
    diff_C1_actual_varied = ...
        (Q_actual_varied(iter,3)- C1_actual )/abs( C1_actual )*100;
    diff_C2_actual_varied = ...
        (Q_actual_varied(iter,4)-C2_actual)/abs( C2_actual )*100;

    fprintf('    Variation Component Set %-2.u: \n', iter);
    fprintf('        R1 = %+11.4f Ohms,    %% diff = %+8.4f (%%).\n', ...
            Q_actual_varied(iter,1), diff_R1_actual_varied);
    fprintf('        R2 = %+11.4f Ohms,    %% diff = %+8.4f (%%).\n', ...
            Q_actual_varied(iter,2), diff_R2_actual_varied );
    fprintf('        C1 = %+11.4e Farads,  %% diff = %+8.4f (%%).\n', ...
            Q_actual_varied(iter,3), diff_C1_actual_varied);
    fprintf('        C2 = %+11.4e Farads,  %% diff = %+8.4f (%%).\n', ...
            Q_actual_varied(iter,4), diff_C2_actual_varied );
    fprintf('            Varied Design H(%+10.4f) = %+8.4f (dB).\n', ...
            f0, Hw_actual_varied_f0(iter));
    fprintf('            Actual Design H(%+10.4f) = %+8.4f (dB).\n', ...
            f0, Hw_actual_f0);
    fprintf('                %% diff = %+8.4f (%%).\n', ...
            diff_0_actual_varied(iter));
    fprintf('            Varied Design H(%+10.4f) = %+8.4f (dB).\n', ...
            f1 , Hw_actual_varied_f1(iter)  );
    fprintf('            Actual Design H(%+10.4f) = %+8.4f (dB).\n', ...
            f1 , Hw_actual_f1  );
    fprintf('                %% diff = %+8.4f (%%).\n', ...
             diff_1_actual_varied(iter) );
end
