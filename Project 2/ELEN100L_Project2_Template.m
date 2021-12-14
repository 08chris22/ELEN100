%% ELEN 100L (Electric Circuits II): Project 2, Christian Garcia Alexander Luo

clear; clc; clf; cla; close all;
format long; format compact;

%% Setup global variables
%

% These Ideal Design element values are fixed in the circuit.
VG  =  1 ;                     % Generator voltage

Vdd_pos  =  15     ;            % Positive power supply voltage
Vdd_neg  = -15     ;            % Negative power supply voltage

R1_ideal_2 = 5000           ;    % Ohms
R2_ideal_2 = 5000           ;    % Ohms
R3_ideal_2 = 400            ;    % Ohms
R4_ideal_2 = 1000           ;    % Ohms
R5_ideal_2 = 1000           ;    % Ohms
C1_ideal_2 = 0.1e-6           ;    % Farads
C2_ideal_2 = 0.1e-6           ;    % Farads

R1_ideal_6 = 1000           ;    % Ohms
R2_ideal_6 = 1000           ;    % Ohms
R3_ideal_6 = 740.938        ;    % Ohms
R4_ideal_6 = 1000           ;    % Ohms
R5_ideal_6 = 1000           ;    % Ohms
C1_ideal_6 = 0.1e-6           ;    % Farads
C2_ideal_6 = 0.1e-6           ;    % Farads

% Build an array for the R elements.
R_ideal_2 = [R1_ideal_2, R2_ideal_2, R3_ideal_2, R4_ideal_2, R5_ideal_2];
R_ideal_6 = [R1_ideal_6, R2_ideal_6, R3_ideal_6, R4_ideal_6, R5_ideal_6];

% Build an array for the C elements.
C_ideal_2 = [ (0),           (0),           (0), (0), (0)         ; ...
              (0), -(C1_ideal_2),           (0), (0), (C1_ideal_2); ...
              (0),           (0), -(C2_ideal_2), (0), (0)         ; ...
              (0),           (0),           (0), (0), (0)         ; ...
              (0),           (0),           (0), (0), (0)         ];

C_ideal_6 = [ (0),           (0),           (0), (0), (0)         ; ...
              (0), -(C1_ideal_6),           (0), (0), (C1_ideal_6); ...
              (0),           (0), -(C2_ideal_6), (0), (0)         ; ...
              (0),           (0),           (0), (0), (0)         ; ...
              (0),           (0),           (0), (0), (0)         ];


% Build an array for the source elements.
B = [ VG; 0; 0; 0; 0 ];

% Build an array for the the time vector.
time2 = [0, 3*10^(-3)];
time6 = [0, 1*10^(-3)];

% Build an array for the the initial conditions.
x0 = [0; 0; 0; 0; 0];     % Assume everything is zero to start

% These values are used for plotting purposes.
fignum = 1;

plot_left_2   = 0;   plot_right_2 = time2(2);   % x-axis range (seconds)
plot_bottom_2 = 0;   plot_top_2   = VG+0.6;     % y-axis range (volts)

plot_left_6   = 0;   plot_right_6 = time6(2);   % x-axis range (seconds)
plot_bottom_6 = 0;   plot_top_6   = VG+0.2;     % y-axis range (volts)


%% Problem 2

% Display the component values for the Ideal design.
%

disp(' ');
disp('The Ideal Design component values are:');
fprintf('    R1 = %+11.4f Ohms.\n', R1_ideal_2   );
fprintf('    R2 = %+11.4f Ohms.\n', R2_ideal_2  );
fprintf('    R3 = %+11.4f Ohms.\n', R3_ideal_2  );
fprintf('    R4 = %+11.4f Ohms.\n', R4_ideal_2  );
fprintf('    R5 = %+11.4f Ohms.\n', R5_ideal_2  );
fprintf('    C1 = %+11.4e Farads.\n', C1_ideal_2 );
fprintf('    C2 = %+11.4e Farads.\n', C2_ideal_2 );

%%
% Calculate the MATLAB transient response for the Ideal design.
%

% Update the resistor variables used in the proj2E100_transient function
% before calling the ode23t solver.
R1_circuit =  R1_ideal_2;
R2_circuit =  R2_ideal_2 ;
R3_circuit =  R3_ideal_2 ;
R4_circuit =  R4_ideal_2 ;
R5_circuit =  R5_ideal_2 ;

options  = odeset('mass', C_ideal_2 , 'RelTol', 0.1e-9);
[t2, x2] = ode23t( @proj2E100_transient , [0 5] , x0 , options );

% Capture peak overshoot and undershoot voltages with indexes.
[v5_pk_overshoot_ideal_2 , ...
    v5_pk_overshoot_ideal_index_2 ] = max( x2(:,5) );
[v5_pk_undershoot_ideal_2, ...
    v5_pk_undershoot_ideal_index_2] = ...
    min( x2(v5_pk_overshoot_ideal_index_2 + 1:size(t2),5) );
v5_pk_undershoot_ideal_index_2 = ...
    v5_pk_undershoot_ideal_index_2 + v5_pk_overshoot_ideal_index_2;

% Capture peak overshoot and undershoot time stamps at peak indexes.
t2_pk_overshoot_ideal_2  = t2( v5_pk_overshoot_ideal_index_2);
t2_pk_undershoot_ideal_2 = t2(v5_pk_undershoot_ideal_index_2);

%%
% Generate the MATLAB plot for the transient response.

fignum = fignum+1; figObj = figure(fignum);  % Establish a figure number
set(fignum, 'Name', ...
    ['Prob 2: Transient Response Ideal Design']); % Name the figure

Tr_ideal_2_Plot = plot(t2, x2);                  % Generate plot
grid on;                                     % Turn grid on
xlabel('Time (seconds)');                    % Label the x-axis
ylabel('Amplitude (volts)');                 % Label the y-axis
axis([plot_left_2, plot_right_2, ...
      plot_bottom_2, plot_top_2]);           % Bound plot
title(['Figure ',num2str(fignum,'%-2.u'),...
       ': Problem 2 Transient Response Ideal Design']);
legend('v_1(t)', 'v_2(t)', 'v_3(t)', 'v_4(t)', 'v_5(t)', ...
    'Location', 'NorthEast');

% Add annotation to the plot.
strmax = ['\leftarrow ', num2str(t2_pk_overshoot_ideal_2 * 1e6), ...
          '\mus, ', num2str(v5_pk_overshoot_ideal_2),'V'];
text(t2_pk_overshoot_ideal_2, ...
     v5_pk_overshoot_ideal_2, ...
     strmax, 'HorizontalAlignment', 'left', 'FontWeight', 'bold');

strmin = ['\leftarrow ', num2str(t2_pk_undershoot_ideal_2 * 1e6), ...
          '\mus, ', num2str(v5_pk_undershoot_ideal_2),'V'];
text(t2_pk_undershoot_ideal_2, ...
     v5_pk_undershoot_ideal_2, ...
     strmin, 'HorizontalAlignment', 'left', 'FontWeight', 'bold');
 
%%
% Display the MATLAB peak overshoot and undershoot values for the Ideal
% design.
%

disp(' ');
disp('The MATLAB peak overshoot and undershoot values are:');
fprintf('    V5 p.o. = %+11.4f Volts.\n',   v5_pk_overshoot_ideal_2);
fprintf('    V5 p.u. = %+11.4f Volts.\n',   v5_pk_undershoot_ideal_2 );
fprintf('     t p.o. = %+11.4e seconds.\n', t2_pk_overshoot_ideal_2 );
fprintf('     t p.u. = %+11.4e seconds.\n', t2_pk_undershoot_ideal_2 );

%% Problem 3
%This part requires screenshots for:
% LTSpice schematic 
% The LTSpice voltage source setup (this is the gray settings box)
% The LT model transient analysis simulation setup(simulation settings)
% The LTSpice simulation result 

% Capture peak overshoot and undershoot voltages from the LTSpice plot.
ltspice_v5_pk_overshoot_ideal_2  =  1.527 ;
ltspice_v5_pk_undershoot_ideal_2 =  0.722 ;

% Capture peak overshoot and undershoot time stamps from the LTSpice plot.
ltspice_t2_pk_overshoot_ideal_2  =  0.000317848 ;
ltspice_t2_pk_undershoot_ideal_2 =  0.00063569 ;

%%
% Display the LTSpice peak overshoot and undershoot values for the Ideal
% design.
%

disp(' ');
disp('The LTSpice peak overshoot and undershoot values are:');
fprintf('    V5 p.o. = %+11.4f Volts.\n', ...
    ltspice_v5_pk_overshoot_ideal_2);
fprintf('    V5 p.u. = %+11.4f Volts.\n', ...
    ltspice_v5_pk_undershoot_ideal_2);
fprintf('     t p.o. = %+11.4e seconds.\n', ...
    ltspice_t2_pk_overshoot_ideal_2);
fprintf('     t p.u. = %+11.4e seconds.\n', ...
    ltspice_t2_pk_undershoot_ideal_2);

%%
% Calculate the percent difference at the peak overshoot and undershoot
% values between MATLAB and LTSpice Ideal Designs.
%

diff_ideal_v5_pk_overshoot_2  = ...
    (ltspice_v5_pk_overshoot_ideal_2 - v5_pk_overshoot_ideal_2) ...
    /abs(v5_pk_overshoot_ideal_2)*100;
diff_ideal_v5_pk_undershoot_2 = (ltspice_v5_pk_undershoot_ideal_2 - v5_pk_undershoot_ideal_2) ...
    /abs(v5_pk_undershoot_ideal_2)*100 ;

diff_ideal_t2_pk_overshoot_2  = (ltspice_t2_pk_overshoot_ideal_2 - t2_pk_overshoot_ideal_2) ...
    /abs(t2_pk_overshoot_ideal_2)*100 ;
diff_ideal_t2_pk_undershoot_2 = (ltspice_t2_pk_undershoot_ideal_2 - t2_pk_undershoot_ideal_2) ...
    /abs(t2_pk_undershoot_ideal_2)*100 ;

disp(' ');
disp('The % difference between MATLAB and LTSpice at the peaks:');
fprintf('    MATLAB  V5 p.o. = %+11.4f Volts.\n', ...
         v5_pk_overshoot_ideal_2);
fprintf('    LTSpice V5 p.o. = %+11.4f Volts.\n', ...
         ltspice_v5_pk_overshoot_ideal_2);
fprintf('        %% diff = %+8.4f (%%).\n', ...
    diff_ideal_v5_pk_overshoot_2);

fprintf('    MATLAB  V5 p.u. = %+11.4f Volts.\n', v5_pk_undershoot_ideal_2 );
fprintf('    LTSpice V5 p.u. = %+11.4f Volts.\n', ltspice_v5_pk_undershoot_ideal_2 );
fprintf('        %% diff = %+8.4f (%%).\n', diff_ideal_v5_pk_undershoot_2 );

fprintf('    MATLAB   t p.o. = %+11.4e seconds.\n', t2_pk_overshoot_ideal_2 );
fprintf('    LTSpice  t p.o. = %+11.4e seconds.\n', ltspice_t2_pk_overshoot_ideal_2 );
fprintf('        %% diff = %+8.4f (%%).\n', diff_ideal_t2_pk_overshoot_2 );
fprintf('    MATLAB   t p.u. = %+11.4e seconds.\n', t2_pk_undershoot_ideal_2 );
fprintf('    LTSpice  t p.u. = %+11.4e seconds.\n', ltspice_t2_pk_undershoot_ideal_2 );
fprintf('        %% diff = %+8.4f (%%).\n', diff_ideal_t2_pk_undershoot_2 );
     
%% Problem 6
% Display the component values for the Ideal design.
%

disp(' ');
disp('The Ideal Design component values are:');
fprintf('    R1 = %+11.4f Ohms.\n',  R1_ideal_6 );
fprintf('    R2 = %+11.4f Ohms.\n',  R2_ideal_6 );
fprintf('    R3 = %+11.4f Ohms.\n',  R3_ideal_6 );
fprintf('    R4 = %+11.4f Ohms.\n',  R4_ideal_6 );
fprintf('    R5 = %+11.4f Ohms.\n',  R5_ideal_6 );
fprintf('    C1 = %+11.4e Farads.\n', C1_ideal_6 );
fprintf('    C2 = %+11.4e Farads.\n', C2_ideal_6 );

%%
% Calculate the MATLAB transient response for the Ideal design.
%

% Update the resistor variables used in the proj2E100_transient function
% before calling the ode23t solver.
R1_circuit =  R1_ideal_6;
R2_circuit =  R2_ideal_6 ;
R3_circuit =  R3_ideal_6 ;
R4_circuit =  R4_ideal_6 ;
R5_circuit =  R5_ideal_6 ;

%ODE solution, refer to the example in Problem 2
options  = odeset('mass', C_ideal_6  , 'RelTol', 0.1e-9);  
[t6, x6] = ode23t( @proj2E100_transient, [0 5] , x0 , options );

% Capture peak overshoot and undershoot voltages with indexes.
[v5_pk_overshoot_ideal_6 , ...
    v5_pk_overshoot_ideal_index_6 ] = max( x6(:,5) );
[v5_pk_undershoot_ideal_6, ...
    v5_pk_undershoot_ideal_index_6] = ...
    min( x6(v5_pk_overshoot_ideal_index_6 + 1:size(t6),5) );
v5_pk_undershoot_ideal_index_6 = ...
    v5_pk_undershoot_ideal_index_6 + v5_pk_overshoot_ideal_index_6;

% Capture peak overshoot and undershoot time stamps at peak indexes.
t6_pk_overshoot_ideal_6  = t6( v5_pk_overshoot_ideal_index_6);
t6_pk_undershoot_ideal_6 = t6(v5_pk_undershoot_ideal_index_6);

%%
% Generate the MATLAB plot for the transient response.

fignum = fignum+1; figObj = figure(fignum);  % Establish a figure number
set(fignum, 'Name', ...
    ['Prob 6: Transient Response Ideal Design']); % Name the figure

Tr_ideal_6_Plot = plot( t6 , x6 );                % Generate plot
grid on;                                     % Turn grid on
xlabel('Time (seconds)');                    % Label the x-axis
ylabel('Amplitude (volts)');                 % Label the y-axis
axis([plot_left_6, plot_right_6, ...
      plot_bottom_6, plot_top_6]);           % Bound plot
title(['Figure ',num2str(fignum,'%-2.u'),...
       ': Problem 6 Transient Response Ideal Design']);
legend('v_1(t)', 'v_2(t)', 'v_3(t)', 'v_4(t)', 'v_5(t)', ...
    'Location', 'NorthEast');

% Add annotation to the plot.
strmax = ['\leftarrow ', num2str(t6_pk_overshoot_ideal_6 * 1e6), ...
          '\mus, ', num2str(v5_pk_overshoot_ideal_6),'V'];
text(t6_pk_overshoot_ideal_6, ...
     v5_pk_overshoot_ideal_6, ...
     strmax, 'HorizontalAlignment', 'left', 'FontWeight', 'bold');

strmin = [num2str(t6_pk_undershoot_ideal_6 * 1e6), ...
          '\mus, ', num2str(v5_pk_undershoot_ideal_6),'V \rightarrow'];
text(t6_pk_undershoot_ideal_6, ...
     v5_pk_undershoot_ideal_6, ...
     strmin, 'HorizontalAlignment', 'right', 'FontWeight', 'bold');
 
%%
% Display the MATLAB peak overshoot and undershoot values for the Ideal
% design.
%

disp(' ');
disp('The MATLAB peak overshoot and undershoot values are:');
fprintf('    V5 p.o. = %+11.4f Volts.\n',   v5_pk_overshoot_ideal_6);
fprintf('    V5 p.u. = %+11.4f Volts.\n',   v5_pk_undershoot_ideal_6);
fprintf('     t p.o. = %+11.4e seconds.\n', t6_pk_overshoot_ideal_6);
fprintf('     t p.u. = %+11.4e seconds.\n', t6_pk_undershoot_ideal_6);

%%
%This part requires screenshots of:
% The LTSpice schematic for the circuit 
% The LTSpice voltage source setup for the circuit 
% The LTSpice model transient analysis simulation setup %%
% The LTSpice simulation result is shown below.

% Capture peak overshoot and undershoot voltages from the plot.
ltspice_v5_pk_overshoot_ideal_6  = 1.09 ;
ltspice_v5_pk_undershoot_ideal_6 = 0.9918 ;

% Capture peak overshoot and undershoot time stamps from the plot.
ltspice_t6_pk_overshoot_ideal_6  = 0.00023978 ;
ltspice_t6_pk_undershoot_ideal_6 = 0.000479554 ;

%%
% Display the LTSpice peak overshoot and undershoot values for the Ideal
% design.
%

disp(' ');
disp('The LTSpice peak overshoot and undershoot values are:');
fprintf('    V5 p.o. = %+11.4f Volts.\n', v5_pk_overshoot_ideal_6 );
fprintf('    V5 p.u. = %+11.4f Volts.\n',  v5_pk_undershoot_ideal_6);
fprintf('     t p.o. = %+11.4e seconds.\n', t6_pk_overshoot_ideal_6 );
fprintf('     t p.u. = %+11.4e seconds.\n', t6_pk_undershoot_ideal_6 );


%%
% Calculate the percent difference at the peak overshoot and undershoot
% values between MATLAB and LTSpice Ideal Designs.
%

diff_ideal_v5_pk_overshoot_6  = ...
    (ltspice_v5_pk_overshoot_ideal_6 - v5_pk_overshoot_ideal_6) ...
    /abs(v5_pk_overshoot_ideal_6)*100;
diff_ideal_v5_pk_undershoot_6 = ...
    (ltspice_v5_pk_undershoot_ideal_6 - v5_pk_undershoot_ideal_6) ...
    /abs(v5_pk_undershoot_ideal_6)*100;

diff_ideal_t6_pk_overshoot_6  = ...
    (ltspice_t6_pk_overshoot_ideal_6 - t6_pk_overshoot_ideal_6) ...
    /abs(t6_pk_overshoot_ideal_6)*100;
diff_ideal_t6_pk_undershoot_6 = ...
    (ltspice_t6_pk_undershoot_ideal_6 - t6_pk_undershoot_ideal_6) ...
    /abs(t6_pk_undershoot_ideal_6)*100;

disp(' ');
disp('The % difference between MATLAB and LTSpice at the peaks:');
fprintf('    MATLAB  V5 p.o. = %+11.4f Volts.\n', ...
         v5_pk_overshoot_ideal_6);
fprintf('    LTSpice V5 p.o. = %+11.4f Volts.\n', ...
         ltspice_v5_pk_overshoot_ideal_6);
fprintf('        %% diff = %+8.4f (%%).\n', ...
         diff_ideal_v5_pk_overshoot_6);
fprintf('    MATLAB  V5 p.u. = %+11.4f Volts.\n', ...
         v5_pk_undershoot_ideal_6);
fprintf('    LTSpice V5 p.u. = %+11.4f Volts.\n', ...
         ltspice_v5_pk_undershoot_ideal_6);
fprintf('        %% diff = %+8.4f (%%).\n', ...
         diff_ideal_v5_pk_undershoot_6);

fprintf('    MATLAB   t p.o. = %+11.4e seconds.\n', ...
         t6_pk_overshoot_ideal_6);
fprintf('    LTSpice  t p.o. = %+11.4e seconds.\n', ...
         ltspice_t6_pk_overshoot_ideal_6);
fprintf('        %% diff = %+8.4f (%%).\n', ...
         diff_ideal_t6_pk_overshoot_6);
fprintf('    MATLAB   t p.u. = %+11.4e seconds.\n', ...
         t6_pk_undershoot_ideal_6);
fprintf('    LTSpice  t p.u. = %+11.4e seconds.\n', ...
         ltspice_t6_pk_undershoot_ideal_6);
fprintf('        %% diff = %+8.4f (%%).\n', ...
         diff_ideal_t6_pk_undershoot_6);

%% Program execution complete
%

disp(' ');
disp('Program execution complete....');

%% MATLAB code listing
%