function [Ts,Tr] = findTsTr(y, t,Nstart, Nend)
    
    yt = y - y(Nstart);
    ss = y(Nend-1) - y(Nstart);
    yt = abs(yt);
    ss = abs(ss);

    Tr = t(Nstart + find(yt(Nstart:Nend) >= ss*0.67, 1)) - t(Nstart);
    yf = flip(yt(Nstart:Nend));
    k = length(yf);
    Ts = t(Nstart + k - find(abs(yf - ss) >= 0.025*ss, 1)) - t(Nstart) ;

end

