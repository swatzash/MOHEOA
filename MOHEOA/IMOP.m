

function fobj = IMOP(name)

    switch name
        case 'IMOP1'
            fobj = @IMOP1;
        case 'IMOP2'
            fobj = @IMOP2; 
        case 'IMOP3'
            fobj = @IMOP3;  
        case 'IMOP4'
            fobj = @IMOP4;
        case 'IMOP5'
            fobj = @IMOP5;
        case 'IMOP6'
            fobj = @IMOP6; 
        case 'IMOP7'
            fobj = @IMOP7;
        case 'IMOP8'
            fobj = @IMOP8; 
    end
end

%% IMOP1
function y = IMOP1(x)
    a1 = 0.05;  % Parameter a1
    K  = 5;     % Parameter K
    
    z1 = mean(x(:,1:K),2).^a1;
    g  = sum((x(:,K+1:end)-0.5).^2,2);
    y(:,1) = g + cos(z1*pi/2).^8;
    y(:,2) = g + sin(z1*pi/2).^8;
    
end

%% IMOP2
function y = IMOP2(x)
    a1 = 0.05;  % Parameter a1
    K  = 5;     % Parameter K

    z1 = mean(x(:,1:K),2).^a1;
    g  = sum((x(:,K+1:end)-0.5).^2,2);
    y(:,1) = g + cos(z1*pi/2).^0.5;
    y(:,2) = g + sin(z1*pi/2).^0.5;
end



%% IMOP3
function y = IMOP3(x)
   
    a1 = 0.05;  % Parameter a1
    K  = 5;     % Parameter K
    z1 = mean(x(:,1:K),2).^a1;
    g  = sum((x(:,K+1:end)-0.5).^2,2);
    y(:,1) = g + (1+cos(z1*pi*10)/5-z1);
    y(:,2) = g + z1;
end

%% IMOP4
function y = IMOP4(x)
    a1 = 0.05;  % Parameter a1
    K  = 5;     % Parameter K 
    z1 = mean(x(:,1:K),2).^a1;
    g  = sum((x(:,K+1:end)-0.5).^2,2);
    y(:,1) = (1+g).*(z1);
    y(:,2) = (1+g).*(z1+sin(10*pi*z1)/10);
    y(:,3) = (1+g).*(1-z1);

end


%% IMOP5
function y = IMOP5(x)
   a1 = 0.05;  % Parameter a1
    a2 = 10;    % Parameter a2
    K  = 5;     % Parameter K  
    z1 = mean(x(:,1:2:K),2).^a1;
    z2 = mean(x(:,2:2:K),2).^a2;
    g  = sum((x(:,K+1:end)-0.5).^2,2);
    y(:,1) = 0.4*cos(pi*ceil(z1*8)/4) + 0.1*z2.*cos(16*pi*z1);
    y(:,2) = 0.4*sin(pi*ceil(z1*8)/4) + 0.1*z2.*sin(16*pi*z1);
    y(:,3) = 0.5 - sum(y(:,1:2),2);
    y = y + repmat(g,1,3);
end


%% IMOP6
function y = IMOP6(x)
   a1 = 0.05;  % Parameter a1
    a2 = 10;    % Parameter a2
    K  = 5;     % Parameter K   

    z1 = mean(x(:,1:2:K),2).^a1;
    z2 = mean(x(:,2:2:K),2).^a2;
    g  = sum((x(:,K+1:end)-0.5).^2,2);
    r  = max(0,min(sin(3*pi*z1).^2,sin(3*pi*z2).^2)-0.05);
    y(:,1) = (1+g).*z1 + ceil(r);
    y(:,2) = (1+g).*z2 + ceil(r);
    y(:,3) = (0.5+g).*(2-z1-z2) + ceil(r);
end



%% IMOP7
function y = IMOP7(x)
    a1 = 0.05;  % Parameter a1
    a2 = 10;    % Parameter a2
    K  = 5;     % Parameter K     
    z1 = mean(x(:,1:2:K),2).^a1;
    z2 = mean(x(:,2:2:K),2).^a2;
    g  = sum((x(:,K+1:end)-0.5).^2,2);
    y(:,1) = (1+g).*cos(z1*pi/2).*cos(z2*pi/2);
    y(:,2) = (1+g).*cos(z1*pi/2).*sin(z2*pi/2);
    y(:,3) = (1+g).*sin(z1*pi/2);
    r = min(min(abs(y(:,1)-y(:,2)),abs(y(:,2)-y(:,3))),abs(y(:,3)-y(:,1)));
    y = y + repmat(10*max(0,r-0.1),1,3);

end

%% IMOP8
function y = IMOP8(x)
    a1 = 0.05;  % Parameter a1
    a2 = 10;    % Parameter a2
    K  = 5;     % Parameter K     

    z1 = mean(x(:,1:2:K),2).^a1;
    z2 = mean(x(:,2:2:K),2).^a2;
    g  = sum((x(:,K+1:end)-0.5).^2,2);
    y(:,1) = z1;
    y(:,2) = z2;
    y(:,3) = (1+g).*(3-sum(y(:,1:2)./(1+repmat(g,1,2)).*(1+sin(19*pi.*y(:,1:2))),2));
end