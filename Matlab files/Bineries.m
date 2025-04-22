function BinOut = Bineries(traceTime, traceVec, VectLimits)
    BinOut = zeros(1,7);
    BinOutAv = zeros(1,7);
    
    if isnan(mean(traceVec(traceTime<=VectLimits(1))))
        AvCntrl = 0;
    else
        AvCntrl = mean(traceVec(traceTime<=VectLimits(1)));
    end;
    
    BinOutAv(1) = AvCntrl;
    
    if AvCntrl~=0
      BinOut(1) = 1;  
    end;
    
    for i=1:length(VectLimits)-1
        if isnan(mean(traceVec((traceTime>=VectLimits(i))&(traceTime<=VectLimits(i+1)))))
            BinOutAv(i+1) = 0;
        else
            BinOutAv(i+1) = mean(traceVec((traceTime>=VectLimits(i))&(traceTime<=VectLimits(i+1))));
        end;
        
        if (BinOutAv(i+1)>=AvCntrl*0.95)&(BinOutAv(i+1)<=AvCntrl*1.05)
            BinOut(i+1) = 0;
        elseif (BinOutAv(i+1)>=AvCntrl*1.05)
            BinOut(i+1) = 1;
        else
            BinOut(i+1) = -1;
        end;
        
    end;
    
end
