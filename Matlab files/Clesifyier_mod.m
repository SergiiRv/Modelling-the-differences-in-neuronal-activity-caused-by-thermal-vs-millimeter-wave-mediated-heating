function ClassOut = Clesifyier_mod(Bx)
    Bx1 = sum(Bx);
    Bx2 = sum(abs(Bx));
    if (Bx1 == 0)&(Bx2 == 0)
        ClassOut = 0; %no AP
    elseif (Bx1 == 1)&(Bx2 == 1)
        ClassOut = -1; % no change
    elseif (Bx1 > 0)&(Bx2 == Bx1)
        ClassOut = 2; % monotonic increase
    elseif  (Bx1 < 0)&((Bx2-2) == (-Bx1))
        ClassOut = -2; % monotonic decrease
    elseif (Bx1 < Bx2)
        ClassOut = 1; % bi-phasic
    else
        ClassOut = -3; % non identified
    end;
end