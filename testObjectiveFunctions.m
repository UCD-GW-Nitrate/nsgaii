function f = testObjectiveFunctions(x, fun)

x = x(:);
n = length(x);
switch fun
    case 'SCH'
        f(1) = x.^2;
        f(2) = (x - 2).^2;
        
    case 'FON'
        a = 1/sqrt(3);
        f(1) = 1 - exp( -sum((x - a).^2) );
        f(2) = 1 - exp( -sum((x + a).^2) );
        
    case 'POL'
        A1 = 0.5*sin(1) - 2*cos(1) + sin(2) - 1.5*cos(2);
        A2 = 1.5*sin(1) - cos(1) + 2*sin(2) - 0.5*cos(2);
        B1 = 0.5*sin(x(1)) - 2*cos(x(1)) + sin(x(2)) - 1.5*cos(x(2));
        B2 = 1.5*sin(x(1)) - cos(x(1)) + 2*sin(x(2)) - 0.5*cos(x(2));
        f(1) = 1 + (A1 - B1)^2 + (A2 - B2)^2;
        f(2) = (x(1) + 3)^2 + (x(2) + 1)^2;
        
    case 'KUR'
        f(1) = sum(-10*exp(-0.2*sqrt(x(1:end-1).^2 + x(2:end).^2)));
        f(2) = sum(abs(x).^0.8 + 5*sin(x.^3));
        
    case 'ZDT1'
        g = 1 + 9*sum(x(2:end))/(n-1);
        f(1) = x(1);
        f(2) = g*(1 - sqrt(x(1)/g));
        
    case 'ZDT2'
        g = 1 + 9*sum(9*x(2:end))/(n-1);
        f(1) = x(1);
        f(2) = g*(1 - (x(1)/g)^2);
        
    case 'ZDT3'
        g = 1 + 9*sum(x(2:end))/(n-1);
        f(1) = x(1);
        f(2) = g*( 1 - sqrt(x(1)/g) - (x(1)/g)*sin(10*pi*x(1)) );
        
    case 'ZDT4'
        g = 1 + 10*(n-1) + sum(x(2:end).^2 - 10*cos(4*pi*x(2:end)));
        f(1) = x(1);
        f(2) = g*(1 - sqrt(x(1)/g));
        
    case 'ZDT6'
        g = 1 + 9*(sum(x(2:end))/(n-1)).^0.25;
        f(1) = 1 - exp(-4*x(1))*sin(6*pi*x(1))^6;
        f(2) = g*(1 - (f(1)/g)^2);
end