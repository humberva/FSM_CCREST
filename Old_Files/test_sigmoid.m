clear; clc;
x = -5:0.001:5;
alp = 1;
plot(x,zeros(length(x),1)+0.5, 'Color', 'k');hold all;
plot(zeros(length(0:0.1:1),1),0:0.1:1, 'Color', 'k');hold all;
cont = 1;
for beta = -2:1:2
    sf = 1./(1+exp(-alp.*(x)+beta));
%     subplot(2,1,1);
    hplot(cont) = plot(x,sf, 'DisplayName', ['b = ', num2str(beta)], 'LineWidth', 2);
    set(gca, 'Xlim', [-5 5], 'Xtick', -5:1:5, 'FontSize', 14);
    hxl = xlabel('x', 'FontSize', 16); hyl = ylabel('f(x)', 'FontSize', 16);
    grid on;
%     subplot(2,1,2);
%     hf = (exp(alp.*(x+beta))-exp(-alp.*(x+beta)))./(exp(alp.*(x+beta))+exp(-alp.*(x+beta)));
%     plot(x,hf); hold all;
%     set(gca, 'Xlim', [-5 5]);
    cont = cont + 1;
end

hleg = legend(hplot, 'Location', 'NorthWest');
set(hleg, 'FontName', 'Symbol', 'FontSize', 16);
set(hxl, 'FontName', 'Times', 'FontAngle', 'Italic', 'FontSize', 16);
set(hyl, 'FontName', 'Times', 'FontAngle', 'Italic', 'FontSize', 16);

print('-depsc', 'Sigmoid_Function_beta_sensitivity.eps');
