function  find_MS_MT(G, C)
T = G/(1+C*G);
S = 1/(1+C*G);

[magT,~,wout] = bode(T);
[magS,~,~] = bode(S, wout);

amT(1, :) = magT(1,1, :);
amS(1, :) = magS(1, 1, :);

[MT, inT] = max(amT);
[MS, inS] = max(amS);
pos = round((1 + length(wout))/2);
figure(gcf);

    p1 = semilogx(wout, amT, 'black');
    hold on;
    p2 = semilogx(wout, amS, 'red--');
        xlabel('log($\omega$) [rad/s]');
        ylabel('$A(\omega)$');
        grid on;
    semilogx(wout(inT), MT, 'b*');
    semilogx(wout(inS), MS, 'b*');
        legend([p1, p2], {'T', 'S'});
        text(wout(round((pos/2))), 1.2, {['Ms = ' num2str(round(MS, 3))], ...
        ['Mt = ' num2str(round(MT, 3))]}, ...
        'Interpreter', 'latex');

end


