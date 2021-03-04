function y = pva(th,pva)

if pva=='p',
    y = [ th^3 th^2 th 1 ];
elseif pva=='v',
    y = [ 3*th^2 2*th 1 0 ];
elseif pva=='a',
    y = [ 6*th 2 0  0 ];
end 