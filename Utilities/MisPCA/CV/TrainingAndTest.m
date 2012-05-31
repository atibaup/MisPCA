function [XTr,Sigma,OmegaTr,OmegaTt]=TrainingAndTest(X,Masks,Options)

S=length(X);

OmegaTr=cell(Options.NCV,1);
OmegaTt=cell(Options.NCV,1);
XTr=cell(Options.NCV,1);
Sigma=cell(Options.NCV,1);
for r=1:Options.NCV
    OmegaTr{r}=cell(S,1);
    OmegaTt{r}=cell(S,1);
    XTr{r}=cell(S,1);
    Sigma{r}=cell(S,1);
    for i=1:S
        p=size(X{i},2);
        n=size(X{i},1);
        OmegaTr{r}{i}=(rand(n,p)>=Options.deltaTest);
        OmegaTt{r}{i}=(ones(n,p)-OmegaTr{r}{i});
        XTr{r}{i}=X{i}.*Masks{i}.*OmegaTr{r}{i};
        p=size(XTr{r}{i},2);
        if p>1
            Sigma{r}{i}=cov(XTr{r}{i}');
        else
            Sigma{r}{i}=XTr{r}{i}*XTr{r}{i}';
        end
        
    end
    
end