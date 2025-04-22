function BinOut = BineriesProb(traceTimeZ, traceVec, VectLimits)
    BinOut = zeros(1,length(VectLimits));
    traceTime = traceTimeZ(traceTimeZ~=0);
    VectLimits_plus = zeros(length(VectLimits)+1,1); VectLimits_plus(1) = 0; VectLimits_plus(2:end) = VectLimits;
    VectLimits_delta = VectLimits_plus(1:end-1)' - VectLimits;
    Proportions = (VectLimits_delta./VectLimits_delta(1)).^-1;
    if isnan(mean(traceVec(traceTime<=VectLimits(1))))
        AvCntrl = 0;
    else
        AvCntrl = numel(traceVec(traceTime<=VectLimits(1)));
    end;
    Total = numel(traceVec(traceTime~=0));
        
    if AvCntrl~=0
      BinOut(1) = AvCntrl/Total;  
    end;
    
    for i=1:length(VectLimits)-1
        if isnan(mean(traceVec((traceTime>=VectLimits(i))&(traceTime<=VectLimits(i+1)))))
            BinOut(i+1) = 0;
        else
            BinOut(i+1) = numel(traceVec((traceTime>=VectLimits(i))&(traceTime<=VectLimits(i+1))))/Total;
        end;               
    end;
   BinOut = BinOut.*Proportions;
end
