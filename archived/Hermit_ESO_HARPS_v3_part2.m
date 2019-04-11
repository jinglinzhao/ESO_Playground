fre         = importdata('frequency.out'); 
FAP_val     = importdata('FAP_50.out'); 
power_fre   = zeros((ORDER+1), length(fre));
for order = 0:ORDER
    power_fre(order+1, :) = importdata(['power', num2str(order), '.out']); 
end
vmax = max(max(power_fre(:, x_idx)));
vmin = FAP_val(1);

order = 1:2:21;
power_order = zeros(length(order), sum(x_idx));
for idx = 1:11
    power_order(idx,:) = mean(abs(coeff(2*idx, :))) * power_fre(2*idx, x_idx);
end
power_nor1 = sum(power_order) / sum(mean(abs(coeff(order+1, :)),2));

order = 0:2:20;
power_order = zeros(length(order), sum(x_idx));
for idx = 1:11
    power_order(idx,:) = mean(abs(coeff(2*idx-1, :))) * power_fre(2*idx-1, x_idx);
end
power_nor2 = sum(power_order) / sum(mean(abs(coeff(order+1, :)),2));
power_max = max(max(power_nor1), max(power_nor2)) * 1.1;

x  = 1./ fre;
x_idx = x < 2000 ;


order = 0:2:20;
ax1 = subplot(80,1,1:28);
    [X_fft,Y_fft] = meshgrid(x(x_idx), order);
    h = pcolor(log(X_fft),Y_fft,   power_fre(order+1, x_idx) );
    set(h, 'edgecolor', 'none');
%     colormap(gray)
    h.FaceColor = 'interp';
    colorbar;
    caxis([vmin vmax]);
    ylabel('Even orders');
    title([star(1:2), ' ', star(3:end)])
    if strcmp(star, 'HD103720')
        xticks([log(4.557), log(10), log(100),log(200)])
        xticklabels({'4.557', '10', '100','200'})    
    else
        xticks([log(10), log(100),log(1000)])
        xticklabels({'10', '100','1000'})  
    end        
%     xticks([log(4.9), log(17.9),log(94), log(217.2)])
%     xticklabels({'4.9', '17.9', '94', '217.2'})
    set(gca,'xticklabel',[])
    yticks(order(1:10))
    set(gca,'fontsize',12)

    
ax2 = subplot(80,1,29:40);
    plot(x(x_idx), power_nor2, 'k', 'LineWidth', 1)
    xlim([min(x(x_idx)) max(x(x_idx))])
    ylim([0, power_max])
    ylabel('Power');
%     yticks(power_nor2)
%     T = text(4.9, 0.11, 'P_1', 'FontSize',12, 'HorizontalAlignment', 'center');
%     T = text(17.9, 0.11, 'P_2', 'FontSize',12, 'HorizontalAlignment', 'center');
%     T = text(217.2, 0.11, 'P_3', 'FontSize',12, 'HorizontalAlignment', 'center');
%     T = text(94, 0.11, 'P*', 'FontSize',12, 'HorizontalAlignment', 'center');   
%     ylim([0 0.1])
    set(gca, 'XScale', 'log')
    set(gca,'fontsize',12)
pos1 = get(ax2, 'Position') % gives the position of current sub-plot
new_pos1 = pos1 +[0.00 0 -0.079 0]
set(ax2, 'Position', new_pos1) % set new po


order = 1:2:21;
ax3 = subplot(80,1,41:52);
    plot(x(x_idx), power_nor1, 'k', 'LineWidth', 1)
    xlim([min(x(x_idx)) max(x(x_idx))])
    ylim([0, power_max])
    ylabel('Power');
%     T = text(4.9, 0.11, 'P_1', 'FontSize',12, 'HorizontalAlignment', 'center');
%     T = text(17.9, 0.11, 'P_2', 'FontSize',12, 'HorizontalAlignment', 'center');
%     T = text(217.2, 0.11, 'P_3', 'FontSize',12, 'HorizontalAlignment', 'center');
%     T = text(94, 0.11, 'P*', 'FontSize',12, 'HorizontalAlignment', 'center');   
%     ylim([0 0.1])
    set(gca,'xticklabel',[])
    set(gca, 'XScale', 'log')
    set(gca,'fontsize',12)
pos1 = get(ax3, 'Position') % gives the position of current sub-plot
new_pos1 = pos1 +[0.00 0 -0.079 0]
set(ax3, 'Position', new_pos1) % set new po
set(gca, 'YDir','reverse')



ax4 = subplot(80,1,53:80);
    [X_fft,Y_fft] = meshgrid(x(x_idx), order);
    h = pcolor(log(X_fft),Y_fft,   power_fre(order+1, x_idx) );
    set(h, 'edgecolor', 'none');
    h.FaceColor = 'interp';
    colorbar;
    caxis([vmin vmax]);
    ylabel('Odd orders');
    if strcmp(star, 'HD103720')
        xticks([log(4.557), log(10), log(100),log(200)])
        xticklabels({'4.557', '10', '100','200'})    
    else
        xticks([log(10), log(100),log(1000)])
        xticklabels({'10', '100','1000'})  
    end
%     xticks([log(4.9), log(17.9),log(94), log(217.2)])
%     xticklabels({'4.9', '17.9', '94','217.2'})
    xlabel('Period [days]');
    yticks(order(1:10))
    set(gca,'fontsize',12)
    set(gca, 'YDir','reverse')



set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 20 30]); %x_width=20cm y_width=30cm
saveas(gcf, [star, 'StackedPperiodogram'], 'png')  