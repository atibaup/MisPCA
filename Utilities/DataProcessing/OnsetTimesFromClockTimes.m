function [OnsetTimes]=OnsetTimesFromClockTimes(Masks,TimeVec,OnsetTimesClock)
S=length(OnsetTimesClock);
for i=1:S
    TimeVecSubj=TimeVec(Masks{i}==1);
    [minVal,minIndx]=min(abs(TimeVec-OnsetTimesClock(i)));
    OnsetTimes(i)=minIndx;
end

