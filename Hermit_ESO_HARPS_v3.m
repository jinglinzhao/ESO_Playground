%%%%%%%%%%
% Update %
%%%%%%%%%%
% integrated odd and even orders in one

%%%%%%%%%%%%%%
% Parameters %
%%%%%%%%%%%%%%
% star        = 'Gl628';
% star        = 'HD103720';
% star        = 'Gl358';
% star        = 'Gl479';
% star        = 'Gl581';
% star        = 'Gl674';
% star        = 'Gl176';
star        = 'Gl388';
MJD         = importdata(['../', star, '/MJD.dat']);
RV_HARPS    = importdata(['../', star, '/RV_HARPS.dat']);
info 		= importdata(['../', star, '/info.dat']);
RVC 		= info(1);
RVW 		= info(2);


cd (['../', star, '/3-ccf_fits/'])
file_list   = dir('*.fits');
file_name   = {file_list.name};
N_FILE      = size(file_name, 2);


ORDER           = 21;                                                        % Highest Hermite order 
array_order     = 0:ORDER;
idx_even        = mod(0:ORDER, 2) == 0;
order_even      = array_order(idx_even);
coeff           = zeros((ORDER+1), N_FILE);
coeff_rvc       = zeros((ORDER+1), N_FILE);
RV_gauss        = zeros(N_FILE,1);


grid_size       = 0.1;
v               = (RVC-RVW : grid_size : RVC+RVW+0.1)';
cd ../../code

%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate Coefficient %
%%%%%%%%%%%%%%%%%%%%%%%%%
h = waitbar(0,'Calculating Gauss-Hermite coefficients for all observations...');
dat_list    = dir(['../', star, '/4-ccf_dat/*.dat']);
dat_name   = {dat_list.name};

for n = 1:N_FILE
    i           = n - 1;
    filename    = ['../', star, '/4-ccf_dat/', char(dat_name(n))];
    A           = importdata(filename);
    f           = fit( v, A, 'a*exp(-((x-b)/c)^2)+d', 'StartPoint', [1/RVW, RVC, RVW/2, 0] );
    b           = RV_HARPS(n);  % shift
    RV_gauss(n) = b;    

    parfor order = 0:ORDER
        temp                    = hermite_nor(order, v - RVC) * grid_size;
        coeff(order+1, n)       = sum(A .* temp);  
    end
    waitbar( n / N_FILE )
end
close(h)  

cd (['../', star, '/'])

[pxx,f] = plomb(coeff(0+1, :), MJD - min(MJD), 0.5,  'normalized');
power = zeros((ORDER+1), size(pxx,1));

for n_hermite = 0:ORDER
    data_write = coeff(n_hermite+1, :)';
    save(strcat('Periodogram_h', sprintf('%02d',n_hermite), '.txt'), 'data_write', '-ascii');         
end
       

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot stacked periodogram %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fre         = importdata('frequency.out'); 
FAP_val     = importdata('FAP_50.out'); 
power_fre   = zeros((ORDER+1), length(fre));
for order = 0:ORDER
    power_fre(order+1, :) = importdata(['power', num2str(order), '.out']); 
end
x  = 1./ fre;
x_idx = x > 0 ; % use all the frequencies returned from the python script. 
order = 0:2:20;
vmax = max(max(power_fre(:, x_idx)));

SMOOTH = 0;

for vmin = [0, FAP_val(1)]

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
    
    ax1 = subplot(80,1,1:28);
        if SMOOTH
            [X_fft,Y_fft] = meshgrid(x(x_idx), order);
            h = pcolor(log(X_fft),Y_fft,   power_fre(order+1, x_idx) );
            h.FaceColor = 'interp';
        else
            [X_fft,Y_fft] = meshgrid(x(x_idx), order-1);
            h = pcolor(log(X_fft),Y_fft,   power_fre(order+1, x_idx) );
        end
%         [X_fft,Y_fft] = meshgrid(x(x_idx), order);
%         h = pcolor(log(X_fft),Y_fft,   power_fre(order+1, x_idx) );
        set(h, 'edgecolor', 'none');
    %     colormap(gray)
%         h.FaceColor = 'interp';
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
        if strcmp(star, 'Gl674')
            xline(32.9,'--r');
            xline(4.6938,'--b');
            T = text(32.9, power_max*0.8, 'P*', 'FontSize',12, 'HorizontalAlignment', 'right');
        end        
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
        if strcmp(star, 'Gl674')
            xline(32.9,'--r');
            xline(4.6938,'--b');
            T = text(4.6938, power_max*0.8, 'P_b', 'FontSize',12, 'HorizontalAlignment', 'right');
        end
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
        
        if SMOOTH
            [X_fft,Y_fft] = meshgrid(x(x_idx), order);
            h = pcolor(log(X_fft),Y_fft,   power_fre(order+1, x_idx) );
            h.FaceColor = 'interp';
        else
            [X_fft,Y_fft] = meshgrid(x(x_idx), order-1);
            h = pcolor(log(X_fft),Y_fft,   power_fre(order+1, x_idx) );
        end
        set(h, 'edgecolor', 'none');
        colorbar;
        caxis([vmin vmax]);
        ylabel('Odd orders');
        if strcmp(star, 'HD103720')
            xticks([log(4.557), log(10), log(100),log(1000)])
            xticklabels({'4.557', '10', '100','1000'})    
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
    if vmin == 0
        saveas(gcf, [star, 'StackedPperiodogram'], 'png')
    else
        saveas(gcf, [star, 'StackedPperiodogram', 'color_min_cut'], 'png')
    end
end
