function [t,x] = rk(x0, dt, tf, f, bt, numinputs)
switch max(size(bt))   % Check Butcher tableau
    case 2 %RK1
        [t,x] = rk1(x0, dt, tf, f, bt, numinputs);
    case 3 %RK2
        [t,x] = rk2(x0, dt, tf, f, bt, numinputs);
    case 5 %RK4
        [t,x] = rk4(x0, dt, tf, f, bt, numinputs);
    otherwise
        disp('Not valid input in rk()!')
        % Error message if not valid Butcher tableau
end

    function [tRK1,xRK1] = rk1(x0, dt, tf, f, bt, numinputs)
        K1 = 0;
        a1 = bt(1,2); % Declare Butcher tableau variables
        b1 = bt(2,2); %
        c1 = bt(1,1); %
        nRK1 = floor(tf / dt);
        tRK1 = zeros(nRK1,numinputs);
        xRK1 = zeros(nRK1,numinputs);
        xRK1(1,:) = x0; % Initial values
        
        for j = 2:nRK1+1
            tRK1(j,:) = tRK1(j-1,:) + dt;
            K1 = f( xRK1(j-1,:) );
            xRK1(j,:) = xRK1(j-1,:) + dt * b1 * K1 ;
        end
        
    end
    function [tRK2,xRK2] = rk2(x0, dt, tf, f, bt, numinputs)
        a11 = bt(1,2);  % Declare Butcher tableau variables
        a12 = bt(1,3);  %
        a21 = bt(2,2);  %
        a22 = bt(2,3);  %
        b1 = bt(3,2);   %
        b2 = bt(3,3);   %
        c1 = bt(1,1);   %
        c2 = bt(2,1);   %
        nRK2 = floor(tf / dt);
        tRK2 = zeros(nRK2,numinputs);
        xRK2 = zeros(nRK2,numinputs);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        xRK2(1,:) = x0;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for j = 2:nRK2+1
            tRK2(j,:) = tRK2(j-1,:) + dt;
            K1 = f( xRK2(j-1,:) );
            K2 = f( xRK2(j-1,:)+a21*dt*K1 );
            xRK2(j,:) = xRK2(j-1,:) + dt * b1 * K1 + dt * b2 * K2 ;
        end
    end
    function [tRK4,xRK4] = rk4(x0, dt, tf, f, bt, numinputs)
        a21 = bt(2,2); a31 = bt(3,2); a32 = bt(3,3); a41 = bt(4,2); a42 = bt(4,3); a43 = bt(4,4);
        b1 = bt(5,2); b2 = bt(5,3); b3 = bt(5,4); b4 = bt(5,5); 
        c1 = bt(1,1); c2 = bt(2,1); c3 = bt(3,1); c4 = bt(4,1);
        
        nRK4 = floor(tf / dt);
        tRK4 = zeros(nRK4,numinputs);
        xRK4 = zeros(nRK4,numinputs);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        xRK4(1,:) = x0;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for j = 2:nRK4+1
            tRK4(j,:) = tRK4(j-1,:) + dt;
            K1 = f( xRK4(j-1,:) );
            K2 = f( xRK4(j-1,:)+a21*dt*K1);
            K3 = f( xRK4(j-1,:)+a31*dt*K1+a32*dt*K2);
            K4 = f( xRK4(j-1,:)+a41*dt*K1+a42*dt*K2+a43*dt*K3);
            xRK4(j,:) = xRK4(j-1,:) + dt*b1*K1 + dt*b2*K2 + dt*b3*K3 + dt*b4*K4;
        end
    end
end