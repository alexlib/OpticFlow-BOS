
function PlotImageGraph(dx,dy,iter,iNum)
%Dx and Dy are the desplacements'matrix 
%iter is the current number of the Optic Flow Iteration
%iNum is the total number of the Optic Flow iterations

[yNum, xNum] = size(dx);
[Y, X]   = ndgrid(1:yNum, 1:xNum);

hold on


% ha = axes('units','normalized','position',[0 0 1 1]);
% uistack(ha,'bottom');
%I = imread('SF_001.jpg'); %Validación flujo supersónico
%I = imread('SRP_stable_type1_crop_001.tif'); %flujo SRP
%I = imread('SRP_stable_type2_crop_001.tif'); %flujo SRP
%I = imread('square_rand_ccc_15_001.tif'); %Gas de encendedor
% I = imread('sop_enc_test_001.tif');

%I = imread('Stable_Heated_001.tif');
%I = imread('HighPressure_001.tif');
I = imread('Unstable_001.tif');
hi_adjust = imadjust(I);
hi = imagesc(I);
colormap gray
% set(ha,'handlevisibility','off','visible','off')
set(hi,'alphadata',.9)



sample  = ceil(yNum/45);
% iter = 100;
% iNum = 100;
IndexY  = 1:sample:yNum;
IndexX  = 1:sample:xNum;
scale   = sample*2;
quiver(X(IndexY,IndexX),Y(IndexY,IndexX),...
               scale*dx(IndexY,IndexX),scale*dy(IndexY,IndexX),'color',[1 1 0]);
title(sprintf('Iteration %d of %d.',iter,iNum));
axis ij equal; axis([-10 xNum+10 -10 yNum+10]);


hold off