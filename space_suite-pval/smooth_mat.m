function smoothed_matrix=smooth_mat(matrix,bin_size, method, revert) 
    
    if nargin<4
        revert=0;
    end
    
    not_explored_full = isnan(matrix);
    matrix(not_explored_full) = 0;
    if method=='gauss'
        smoothed_matrix=filter(gausswin(2),1,matrix/max(matrix));
        if sum(isnan(smoothed_matrix))>0
            smoothed_matrix=matrix;
        end
    else   
        kernel_size = [bin_size bin_size];
        occupancy_std = 2;

        [Xgrid,Ygrid]=meshgrid(-kernel_size(1)/2: kernel_size(1)/2, -kernel_size(2)/2:kernel_size(2)/2);
        Rgrid=sqrt((Xgrid.^2+Ygrid.^2));
        kernel = pdf('Normal', Rgrid, 0, occupancy_std);
        kernel = kernel./sum(sum(kernel));
        smoothed_matrix = conv2(matrix, kernel, 'same'); % smoothing
    end
    if revert==1
        smoothed_matrix(not_explored_full) = NaN;
    end
end