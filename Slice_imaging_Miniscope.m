clear all; close all

time_frequency = 10 %Hz 
xval = 1/time_frequency;
x = xval:xval:1800; %last time

F = xlsread("file name", sheet number, matrix);
Window = [100,500]; %Adjust to your conditions

windnum=size(Window,2);
[rownum,colnum]=size(F);

%F is the raw fluorescence signal
%W is the width of the sliding Window (~300 seems good for your data).
%But play with it to see how it changes
%%

for ii = 1:windnum
    W=Window(ii)
    N = length(F);
    for i = 1:N  %Loop through each time sample
        
        dom = i-W/2:i+W/2;
        dom = dom(find(dom>=1 & dom<=N));
        pc = F(dom);  %Get a window of the fluorescence signal around time point i
        mi = prctile(pc,5);
        ma= prctile(pc,15);
        id = find(pc>mi & pc<ma);
        F0(i) = median(pc(id));  %To get F0 at this time point, compute the median based on a subsample of the window
        
    end
    
    subplot(2,windnum,ii)
    plot(x, F)
    hold on
    plot(x, F0,'r')
    legend('F','F0')
    
    subplot(2,windnum,windnum+ii)
    dF=(F-F0')./F0'*100;
    plot(x, dF)
    ylabel('dF/F0(%)')
end
