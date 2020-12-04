function [ y ] = rTeagerCompute( x,VL )
    % calculates Teager energy for the input signal
    % 
    % INPUT: x -> 1xk vector 
    %        VL -> averaging window size, VL=2 for point-wise Teager
    
    size_x=size(x);
    if(rem(VL,2))
        disp('error! VL should be even!');
        exit
    end
    k=length(x);
    start_point=VL/2;
    yy=zeros(size_x);
    for i=1:start_point
        temp_a=x(start_point+1:k-start_point).*conj(x(start_point+1:k-start_point))...
            -x(start_point+1+i:k-start_point+i).*conj(x(start_point+1-i:k-start_point-i));
        yy(start_point+1:k-start_point)=yy(start_point+1:k-start_point)+temp_a;
    end
    yy(1:start_point) = yy(start_point+1);
    yy(k-start_point+1:k) = yy(k-start_point);
    % y=yy(start_point+1:k-start_point);
    y=abs(yy);
    
    % OUTPUT: y -> 1xk vector containg Teager values  
end