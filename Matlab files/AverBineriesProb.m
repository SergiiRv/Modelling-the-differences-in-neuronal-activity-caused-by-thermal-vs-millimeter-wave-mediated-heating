function BinOut = AverBineriesProb(traceTimeZ, traceVec, VectLimits)
    BinOut = zeros(1,length(VectLimits));
    traceTime = traceTimeZ(traceTimeZ~=0);   
    if isnan(mean(traceVec(traceTime<=VectLimits(1))))
        AvCntrl = 0;
    else
        AvCntrl = mean(traceVec(traceTime<=VectLimits(1)));
    end;
    
    if isnan(mean(traceVec((traceTime>=VectLimits(end - 1))&(traceTime<=VectLimits(end)))))
            BinOut_ = 0;
    else
            BinOut_ = mean(traceVec((traceTime>=VectLimits(end - 1))&(traceTime<=VectLimits(end))));
    end; 
    
    if AvCntrl~=0
      BinOut = BinOut_/AvCntrl;
    else
       BinOut = BinOut_/AvCntrl;%5;
    end;
    
end