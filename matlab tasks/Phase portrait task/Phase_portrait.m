function Phase_portrait( func, time,  xlimit, ylimit, xdotsnum, ydotsnum)
%Phase_portrait draws phase portrait of system declared by func
%   for the time  from zeroes to xlimit and ylimit with counts of dots
%   delared by xdotsnum and ydotsnum
        ystep = ylimit/ydotsnum;
        xstep = xlimit/xdotsnum;
      
        
        for y0=ylimit-2:ystep:ylimit
            for x0=xlimit-2:xstep:xlimit
                if (x0 == xlimit-2)||(x0 == xlimit)||(y0 == ylimit)||(y0 == ylimit-2)
                 
                    [T, Y] = ode45(func, time, [x0; y0]);%FUNCTION
                    x = Y(:, 1);
                    y = Y(:, 2);
                   % plot(x(1), y(1), 'go')%starting point
                    plot(x, y, 'b')
                    %plot(x(numel(x)), y(numel(y)), 'r*')%ending point
                
            end
        end

end

