clc
close all
clear
y_ST=[
253.144
993.692
2698.864
2858.6
3613.576
];
y_LO=[
301.356
996.564
2751
3091.692
3476.012
];
x = [1:5];
xlabelnames={
    '(10,3)'
    '(20,3)'
    '(30,3)'
    '(40,3)'
    '(50,3)'
    };
figure
plot(x,y_ST,'-o',x,y_LO,'-.o')
title('Average run time for different instance sizes')
xlabel('Instances size') % x-axis label
ylabel('Average run time')
legend('Strict TW','Loose TW')
set(gca,'xtick',1:7); 
set(gca,'xticklabel',xlabelnames);
set(gca,'XTickLabelRotation',30);
