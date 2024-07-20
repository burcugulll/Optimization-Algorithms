function [obj, tahminEdilenKumeMerkezi] = objfunc_final(pos, veri)


minUzaklik = 0;  

maliyet = 0;  

    for i = 1:size(veri, 1)

        d1 = sqrt(sum((pos(1, 1:7) - veri(i, :)).^2));
        d2 = sqrt(sum((pos(1, 8:14) - veri(i, :)).^2));
        d3 = sqrt(sum((pos(1, 15:21) - veri(i, :)).^2));

        [minUzaklik, indis] = min([d1 d2 d3]);

        maliyet = maliyet+minUzaklik;
        tahminEdilenKumeMerkezi(i,1)= indis;
    end

    obj = maliyet;
end
