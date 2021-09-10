e=zeros(3,1);

calculated=3;
reference=2.5;
threshold=1;

e(1)=abs(calculated-reference)>threshold; % error1 occured
e(2)=0; % error2 not occured
e(3)=0; % error3 not occured
%% test if parallelization can be used
parfor i=1:3
    a=i;
end
%%
fprintf('error_status \n');
disp(e);

if sum(e)>=1

    error('error occured in the test');
end