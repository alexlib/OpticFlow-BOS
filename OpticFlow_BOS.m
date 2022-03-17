clc
clear all
close all
tic

%% Set files location

%%Validation images Supersonic Flow
filePattern = 'C:/Users/Jesica González/Documents/School/UNAM/Titulación/Thesis/CódigoFlujoOptico-BOS/OpticFlow-BOS/ValidationImagesSF/SF_%03d.jpg';


%%Images for optic flow vs cross-correlation tests

%Linear displacement
%filePattern = 'C:/Users/Jesica González/Documents/School/UNAM/Titulación/Thesis/CódigoFlujoOptico-BOS/OpticFlow-BOS/L_3p_2d_%03d.tif';

%Vortex displacements
%filePattern = 'C:/Users/Jesica González/Documents/School/UNAM/Titulación/Thesis/CódigoFlujoOptico-BOS/OpticFlow-BOS/V_3p_2d_%03d.tif';


%%Pattern Tests

%Density of dots
%filePattern = 'C:/Users/Jesica González/Documents/School/UNAM/Titulación/Thesis/CódigoFlujoOptico-BOS/OpticFlow-BOS/Rand_pattern_black_crop_600000_%03d.tif';
%filePattern = 'C:/Users/Jesica González/Documents/School/UNAM/Titulación/Thesis/CódigoFlujoOptico-BOS/OpticFlow-BOS/Rand_pattern_black_crop_60000_%03d.tif';

%Structure of dots
%filePattern = 'C:/Users/Jesica González/Documents/School/UNAM/Titulación/Thesis/CódigoFlujoOptico-BOS/OpticFlow-BOS/square_rand_15_%03d.tif';
%filePattern = 'C:/Users/Jesica González/Documents/School/UNAM/Titulación/Thesis/CódigoFlujoOptico-BOS/OpticFlow-BOS/square_15_%03d.tif';


%%Object of study

%Gas of lighter 15cm
%filePattern = 'C:/Users/Jesica González/Documents/School/UNAM/Titulación/Thesis/CódigoFlujoOptico-BOS/OpticFlow-BOS/square_rand_ccc_15_%03d.tif';

%Blowtorch
%filePattern = 'C:/Users/Jesica González/Documents/School/UNAM/Titulación/Thesis/CódigoFlujoOptico-BOS/OpticFlow-BOS/sop_enc_test_%03d.tif';


%% Read Image Sequence

ImgSeq = readImgSeq(filePattern,0,1);


%% Estimate optic flow for the image sequence

lx=1;
ly=1;

opt.eta = 0.1;
[Matx, Maty] = estimateOpticFlow2D(ImgSeq,opt); %Optic flow displacements' matrixes


%% Scale displacements' data using cross-correlation results

data_PIV_x = load('u_PIV_validation2.mat'); %Supersonic Flow Validation Images
data_PIV_y = load('v_PIV_validation2.mat');

% data_PIV_x = load('u_PIV_square_rand_ccc_15.mat'); %Gas of Lighter
% data_PIV_y = load('v_PIV_square_rand_ccc_15.mat');

% data_PIV_x = load('u_PIV_sop_enc_test.mat'); %Gas of blowtorch
% data_PIV_y = load('v_PIV_sop_enc_test.mat');

data_PIV_x = data_PIV_x.u_original;
data_PIV_y = data_PIV_y.v_original;
data_PIV_x = data_PIV_x{1,1};
data_PIV_y = data_PIV_y{1,1};

%Set last window and step size
window = 16; %For Supersonic Flow, square_rand_crop_15 and sop_enc_test
step = 8;

n = size(data_PIV_y,2); %Number of columns of cross-correlation matrix

u = zeros(size(Matx)); %Generate initial matrixes for scaled displacements (u & v)
v = zeros(size(Maty));

diff_limit_lower = 0.05; %limit for displacements < 1 px
diff_limit_upper = 0.5; %limit for displacements > 1 px

for ii = 1:n %Work with each column from cross-correlation displacements' matrixes
    
    %Define limits of columns subset for the opfic flow matrixes
    left_s = (ii*step) - ((step/2)-1);
    right_s = (ii*step) + (step/2);
    
    for w = left_s:right_s %Work with each column from the optic flow subset
        
        %Asign scaled displacements' magnitudes to each column of u and v
        [u(:,w),v(:,w)] = Scaling(Matx(:,w),Maty(:,w), data_PIV_x(:,ii),data_PIV_y(:,ii), diff_limit_lower, diff_limit_upper,step);
    
    end
    
end


%% BOS Analysis

%Experimental parameters

%Gas of lighter and blowtorch

% esc = 22222.22;             %Inverse value of scale = 1/m
% h = 0.0015;                 %Size of object of study - Lighter
% L = 0.15;                   %Bakground/pattern - object of study distance
% Mag = 0.0858;               %Magnification 
% n_0=1.0002921;              %Refractive Index of sorrounding gas (air)
% rho_0=1.204;                %Density of sorrounding gas (air)
% G = 5.3236e-4;              %Gladstone-Dale constant for butane
% scale= (1/esc);             %Scale -> 1 px = X meters


%Supersonic Flow

esc = 47167;                  %Inverse value of scale = 1/m
h = 0.004;                    %Size of object of study - Nozzle
L = 0.147;                    %Bakground/pattern - object of study distance
Mag = 0.94;                   %Magnification 
n_0=1.0002921;                %Refractive Index of sorrounding gas (air)
rho_0=1.204;                  %Density of sorrounding gas (air)
G = 2.2649e-4;                %Gladstone-Dale constant for air
scale= (1/esc);               %Scale -> 1 px = X meters



Matxx=u; 
Matyy=v; 

%Set values for Finitie Differences routine
m=2; %number of consecutive pixel cells to analyze
a=1; %Distance between pixels
oo=1;

Dm=DifFinFun_2(Matxx,Matyy,m,a)/esc; %Finite Differences matrix

%Set values for Poisson Solver routine
lx=size(Dm,1); 
ly=size(Dm,2);

k= (n_0)/(G*Mag*L*h); %Constant term on Poisson equation
RHS=real((k)*Dm); %

[Density] = Poisson(real(RHS),lx,ly,m,rho_0,oo,a); %Compute density matrix

[IRefraction]= Density*G+1; %Compute density matrix

%Graphs

figure11 = figure;
axes1 = axes('Parent',figure11);
mesh_meters(Density,1,scale)
shading interp
view(0,-90)
colorbar('peer',axes1);

figure22 = figure;
axes2 = axes('Parent',figure22);
mesh_meters(IRefraction,1,scale)
shading interp
view(0,-90)
colorbar('peer',axes1);



toc