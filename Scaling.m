function [data_x, data_y] = Scaling(data_Optic_flow_x,data_Optic_Flow_y, data_CC_x,data_CC_y, diff_limit_lower, diff_limit_upper,step)

% xdiff_limit refres to the max difference between optic flow and cross-correlation (CC) data 
% This function scale the data from each column in the optic flow subset and applies a three point mean to smooth the data
% All data sets need to be of size mx1

m = size(data_CC_x,1); %Number of row in each column

data_x = zeros(size(data_Optic_flow_x,1),1); %Generate initial matrixes for scaled displacements
data_y = zeros(size(data_Optic_Flow_y,1),1);

data_interval = 0.1; %Set data interval
scale_x = General_scale(data_Optic_flow_x, data_CC_x,data_interval); %Get general scale for current column (horizontal displacements' matrix)
scale_y = General_scale(data_Optic_Flow_y, data_CC_y,data_interval); %Get general scale for current column (vertical displacements' matrix)


%% Relate each cross-correlation column data to every optic flow row data 

for k = 1:m %k value in the current cross-correlation column
    
    clear left_val_1;
    clear left_val_2;
    clear right_val_1;
    clear right_val_1;
    
    
    %Define range for optic flow subset of rows 
    
    left_val_1 = (k*step)-((step/2)-1);
    right_val_2 = k*step+(step/2);
    
    right_val_1 = k*step;
    left_val_2 = (k*step)+1;
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%X VALUES%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %Substitute any NaN data with the pixel neighbors' mean value
    if isnan(data_CC_x(k,1))
        if k == 1 
            if isnan(data_CC_x(k+1,1))
                data_CC_x(k,1) = data_CC_x(k+2,1);
            else
                data_CC_x(k,1) = data_CC_x(k+1,1);
            end
        elseif k == m
            data_CC_x(k,1) = data_CC_x(k-1,1); 
        elseif k== m-1
            if isnan(data_CC_x(k+1,1)) && ~isnan(data_CC_x(k-1,1))
                data_CC_x(k,1) = data_CC_x(k-1,1);
            elseif isnan(data_CC_x(k-1,1)) && ~isnan(data_CC_x(k+1,1))
                data_CC_x(k,1) = (data_CC_x(k-2,1)+ data_CC_x(k+1,1))/2;
            else
                data_CC_x(k,1) = (data_CC_x(k-1,1)+ data_CC_x(k+1,1))/2;
            end
        elseif k== m-2
            if isnan(data_CC_x(k+1,1)) && ~isnan(data_CC_x(k-1,1))
                data_CC_x(k,1) = data_CC_x(k-1,1);
            elseif isnan(data_CC_x(k-1,1)) && ~isnan(data_CC_x(k+1,1))
                data_CC_x(k,1) = (data_CC_x(k-2,1)+ data_CC_x(k+1,1))/2;
            else
                data_CC_x(k,1) = (data_CC_x(k-1,1)+ data_CC_x(k+1,1))/2;
            end
        else
            if isnan(data_CC_x(k+1,1)) && ~isnan(data_CC_x(k+2,1))
                data_CC_x(k,1) = (data_CC_x(k-1,1)+ data_CC_x(k+2,1))/2;
                %data_PIV_x(k,1) = data_PIV_x(k-1,1);
            elseif isnan(data_CC_x(k-1,1)) && ~isnan(data_CC_x(k+1,1))
                data_CC_x(k,1) = (data_CC_x(k-2,1)+ data_CC_x(k+1,1))/2;
            elseif ~isnan(data_CC_x(k-1,1)) && isnan(data_CC_x(k+1,1))
                data_CC_x(k,1) = (data_CC_x(k-1,1));
            else
                data_CC_x(k,1) = (data_CC_x(k-1,1)+ data_CC_x(k+1,1))/2;
            end
        end
    end
        
        
    %Evaluate each optic flow row value (i) with current cross-correlation value (k)
    
    for i = left_val_1:right_val_2 %Extra condition considering Horn's set of 8 values
        diff_Horn_PIV_x = abs(data_CC_x(k,1) - data_Optic_flow_x(i,1)); %Difference between optic flow and cross-correlation

        if diff_Horn_PIV_x > diff_limit_lower %Scaling is needed
            n_interval_x = data_Optic_flow_x(i,1)/data_interval;
            data_x(i,1) = data_Optic_flow_x(i,1) + n_interval_x*scale_x*0.4;%Proposed scaling


            %Sum or sustract the diference
            if data_CC_x(k,1) < 0
                data_x(i,1) = data_x(i,1) - diff_Horn_PIV_x;
            else
                data_x(i,1) = data_x(i,1) + diff_Horn_PIV_x;
            end     

        else %No scaling is needed
            data_x(i,1) = data_Optic_flow_x(i,1);
        end 
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%Y VALUES%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %Substitute any NaN data with the pixel neighbors' mean value
    if isnan(data_CC_y(k,1))
        if k == 1 
            if isnan(data_CC_y(k+1,1))
                data_CC_y(k,1) = data_CC_y(k+2,1);
            else
                data_CC_y(k,1) = data_CC_y(k+1,1);
            end
        elseif k == m
            data_CC_y(k,1) = data_CC_y(k-1,1);
           
        elseif k== m-1
            if isnan(data_CC_y(k+1,1)) && ~isnan(data_CC_y(k-1,1))
                data_CC_y(k,1) = data_CC_y(k-1,1);
            elseif isnan(data_CC_y(k-1,1)) && ~isnan(data_CC_y(k+1,1))
                data_CC_y(k,1) = (data_CC_y(k-2,1)+ data_CC_y(k+1,1))/2;
            else
                data_CC_y(k,1) = (data_CC_y(k-1,1)+ data_CC_y(k+1,1))/2;
            end
        elseif k== m-2
            if isnan(data_CC_y(k+1,1)) && ~isnan(data_CC_y(k-1,1))
                data_CC_y(k,1) = data_CC_y(k-1,1);
            elseif isnan(data_CC_y(k-1,1)) && ~isnan(data_CC_y(k+1,1))
                data_CC_y(k,1) = (data_CC_y(k-2,1)+ data_CC_y(k+1,1))/2;
            else
                data_CC_y(k,1) = (data_CC_y(k-1,1)+ data_CC_y(k+1,1))/2;
            end
        else
            if isnan(data_CC_y(k+1,1)) && ~isnan(data_CC_y(k+2,1))
                data_CC_y(k,1) = (data_CC_y(k-1,1)+ data_CC_y(k+2,1))/2;
                %data_PIV_x(k,1) = data_PIV_x(k-1,1);
            elseif isnan(data_CC_y(k-1,1)) && ~isnan(data_CC_y(k+1,1))
                data_CC_y(k,1) = (data_CC_y(k-2,1)+ data_CC_y(k+1,1))/2;
            elseif ~isnan(data_CC_y(k-1,1)) && isnan(data_CC_y(k+1,1))
                data_CC_y(k,1) = (data_CC_y(k-1,1));
            else
                data_CC_y(k,1) = (data_CC_y(k-1,1)+ data_CC_y(k+1,1))/2;
            end
        end
    end
      
    
 
        
    for i = left_val_1:right_val_2
        diff_Horn_PIV_y = abs(data_CC_y(k,1) - data_Optic_Flow_y(i,1)); %Difference between optic flow and cross-correlation 

        if diff_Horn_PIV_y > diff_limit_lower %Scaling is needed
            n_interval_y = data_Optic_Flow_y(i,1)/data_interval;
            data_y(i,1) = data_Optic_Flow_y(i,1) + n_interval_y*scale_y*0.4;%Proposed scaling

            %Sum or sustract the diference
            if data_CC_y(k,1) < 0
                data_y(i,1) = data_y(i,1) - diff_Horn_PIV_y;
            else
                data_y(i,1) = data_y(i,1) + diff_Horn_PIV_y;
            end     

        else %No scaling is needed
            data_y(i,1) = data_Optic_Flow_y(i,1);
        end   
    end
        
end


%% Fix beginning and end of each subset on the column
        
datax_mean = data_x;
datay_mean = data_y;

for k_mean = 1:m-1
    
    left_val = (k_mean*step)-((step/2)-1); 
    right_val = k_mean*step+(step/2);

    %X data matrix
    datax_mean(left_val,1) = (data_x(left_val-1,1) + data_x(left_val,1) + data_x(left_val+1,1))/3;
    datax_mean(right_val,1) = (data_x(right_val-1,1) + data_x(right_val,1) + data_x(right_val+1,1))/3;
    
    %Y data matrix
    datay_mean(left_val,1) = (data_y(left_val-1,1) + data_y(left_val,1) + data_y(left_val+1,1))/3;
    datay_mean(right_val,1) = (data_y(right_val-1,1) + data_y(right_val,1) + data_y(right_val+1,1))/3;
    
end
