function scale = General_scale(data_Optic_Flow, data_CC,data_interval)

%data_Optic_Flow and data_CC must be vectors of mx1


[m,n]=size(data_Optic_Flow);

interval_Horn = max(data_Optic_Flow)-min(data_Optic_Flow);
interval_PIV = max(data_CC)-min(data_CC);


scale = data_interval*(interval_PIV/interval_Horn); 