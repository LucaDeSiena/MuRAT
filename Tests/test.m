e=zeros(3,1);

calculated=3;
reference=2;
threshold=0.5;
e(1)=abs(calculated-reference)>threshold; % error1 occured
e(2)=0; % error2 not occured
e(3)=0; % error3 occured
e(1)=0;

fprintf('test result \n');
disp(e)

if sum(e)>=1
    error('error occured in the test');
end