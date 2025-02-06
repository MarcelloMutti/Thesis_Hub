clear all; clc;

load("VP_res\VP_2014FO.mat")
load("VP_res\VP_2014TO.mat")

tfS=[];
t0S=[];
t0sk=[];
tfsk=[];

% load("VP_res\VP_2014_sample_A.mat")
% tfS=[tfS, sol.tf_ad];
% t0S=[t0S, sol.t0];
% load("VP_res\VP_2014_sample_B.mat")
% tfS=[tfS, sol.tf_ad];
% t0S=[t0S, sol.t0];
% load("VP_res\VP_2014_sample_C.mat")
% tfS=[tfS, sol.tf_ad];
% t0S=[t0S, sol.t0];
% load("VP_res\VP_2014_sample_D.mat")
% tfS=[tfS, sol.tf_ad];
% t0S=[t0S, sol.t0];

%--------------------------------------------------------------------------

cspice_furnsh('kernels\metaker.tm')
fprintf('\nTOTAL kernels number: %d\n', cspice_ktotal('ALL'));

LU=cspice_convrt(1,'AU','KM');              % 1AU [km]
TU=sqrt(LU^3/cspice_bodvrd('Sun','GM',1));  % mu_S=1

if ~isempty(t0S)
    t0S=et2MJD2000(t0S);
    tfS=tfS*TU/86400;
end

TO_sk=EO_prob(strcmp({EO_prob.sts},'skp'));
if ~isempty(TO_sk)
    t0sk=et2MJD2000([TO_sk.t0]);
    tfsk=[TO_sk.tf_ad]*TU/86400;
end

t0 = et2MJD2000([EO_prob.t0]); % Example initial times
tf = [EO_prob.tf_ad]*TU/86400; % Example final times
mp = [EO_prob.mp];  % Example propellant masses

% Create a regular grid
N=500;
t0_grid = linspace(min(t0), max(t0), N);
tf_grid = linspace(min(tf), max(tf), N);
[T0, TF] = meshgrid(t0_grid, tf_grid);

% Interpolate mp values onto the grid
mp_grid = griddata(t0, tf, mp, T0, TF, 'linear'); % Use 'linear', 'cubic', or 'nearest'

% Define the boundary using the boundary function
k = boundary(t0.', tf.', 0.9); % Adjust shrink factor if needed (e.g., 0.8 for tighter fit)

% Create a mask for points inside the boundary
[inBoundary, onBoundary] = inpolygon(T0, TF, t0(k), tf(k)); % Check points inside or on the boundary

% Mask the points outside the boundary
mp_grid(~inBoundary) = NaN; % Set values outside the boundary to NaN

% Set default font size for axes labels, tick labels and legend
fs=15;
set(0,'DefaultAxesFontSize',fs);
set(0,'DefaultLegendFontSize',fs);

% Set the default linewidth
set(0,'DefaultLineLineWidth',2);

% Set the default text and legend interpreter
set(groot,'defaultTextInterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

% Plot the colormap
figure;
contourf(T0, TF, mp_grid, 0:0.2:max([EO_prob.mp]), 'fill', 'on')
c=colorbar;
hold on
contour(T0, TF, mp_grid, [2.8 2.8], 'k--', 'LineWidth', 2)
plot(et2MJD2000([TO_prob.t0]),TU/86400*([TO_prob.tf_ad]),'r','LineWidth',2)
if ~isempty(t0sk)
    plot(t0sk,tfsk,'sk')
end
if ~isempty(t0S)
    plot(t0S,tfS,'ok')
    text(t0S(1)-35,tfS(1)-35,'A','Interpreter','latex','FontSize',fs)
    text(t0S(2)-35,tfS(2)-35,'B','Interpreter','latex','FontSize',fs)
    text(t0S(3)-35,tfS(3)-35,'C','Interpreter','latex','FontSize',fs)
    text(t0S(4)-35,tfS(4)-35,'D','Interpreter','latex','FontSize',fs)
end
hold off
grid on
grid minor
ylim([100,800])
xlabel('$$t_0\, \left[MJD2000\right]$$','Interpreter','latex')
ylabel('$$ToF\, \left[days\right]$$','Interpreter','latex')
c.Label.Interpreter = 'latex';
c.Label.String='$$m_p\, \left[kg\right]$$';
c.Label.FontSize = fs;

% User plot setting removal
set(0,'DefaultAxesFontSize','remove');
set(0,'DefaultLegendFontSize','remove');
set(0,'DefaultLineLineWidth','remove');
set(groot,'defaultTextInterpreter','remove');
set(groot,'defaultLegendInterpreter','remove');

fprintf('Total Solutions: %.0f\n',length(EO_prob))
fprintf('Total ToF skips: %.2f\n',100*(1-length(TO_sk)/length(TO_prob)))
fprintf('Total E2F skips: %.2f\n',100*(1-sum([EO_prob.epsilon]~=0)/length(EO_prob)))