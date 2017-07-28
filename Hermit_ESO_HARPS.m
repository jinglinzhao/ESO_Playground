% Branch from Hermit_coeff_nor.m

% Integrated from Hermit_coeff_NOR_0720_GI479.m

%%%%%%%%%%
% Update %
%%%%%%%%%%


%%%%%%%%%%%%%%
% Parameters %
%%%%%%%%%%%%%%
cd ../HD224817/3-ccf_fits/
file_list   = dir('*.fits');
file_name   = {file_list.name};
N_FILE      = size(file_name, 2);
MJD         = zeros(N_FILE, 1);

for n = 1:N_FILE
    FILE    = file_name{n};
    spec    = fitsread(char(FILE));                                         % spectra of all orders
    H       = fitsinfo(char(FILE));                                         % fits information including header
    header  = H.PrimaryData.Keywords;                                       % get header only
    MJD(n)  = header{30, 2};
end

ORDER           = 5;                                                        % Highest Hermite order 
array_order     = 0:ORDER;
idx_even        = mod(0:ORDER, 2) == 0;
order_even      = array_order(idx_even);
coeff           = zeros((ORDER+1), N_FILE);
coeff_rvc       = zeros((ORDER+1), N_FILE);
RV_gauss        = zeros(N_FILE,1);


grid_size       = 0.1;
v               = (-20 : grid_size : -2)';                                  % km/s
cd ../../code

%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate Coefficient %
%%%%%%%%%%%%%%%%%%%%%%%%%
h = waitbar(0,'Calculating Gauss-Hermite coefficients for all observations...');
dat_list    = dir('../HD224817/4-ccf_dat/*.dat');
dat_name   = {dat_list.name};

for n = 1:N_FILE
    i           = n - 1;
    filename    = ['../HD224817/4-ccf_dat/', char(dat_name(n))];
    A           = importdata(filename);
    f           = fit( v, A, 'a*exp(-((x-b)/c)^2)+d', 'StartPoint', [0.15, -11.28, 2.5, 0] );
    b           = f.b;  % shift
    % plot(v,A, v,f(v)) % test %
    RV_gauss(n) = b;    

    for order = 0:ORDER
        temp                    = hermite_nor(order, v + 11.28) * grid_size;
        temp_rvc                = hermite_nor(order, v - b) * grid_size;
        coeff(order+1, n)       = sum(A .* temp);  
        coeff_rvc(order+1, n)   = sum(A .* temp_rvc); 
    end
    waitbar( n / N_FILE )
end
close(h)  

cd ..

for order = 0:ORDER
    
    [pxx,f] = plomb(coeff(order+1, :), MJD - min(MJD), 0.1, 100, 'normalized');
    [pmax,lmax] = max(pxx);
    f0 = f(lmax);
    disp(['T_planet: ', num2str(1/f0)]);

    [pxx_rvc,f_rvc] = plomb(coeff_rvc(order+1, :), MJD - min(MJD), 0.1, 100, 'normalized');
    [pmax_rvc,lmax_rvc] = max(pxx_rvc);
    f0_rvc = f_rvc(lmax_rvc);
    disp(['T_activity: ', num2str(1/f0_rvc)]);

    [pks,locs] = findpeaks(pxx, f);                 % find all the peaks in (pxx, f)
    [pks_maxs, idx_maxs] = sort(pks, 'descend');    % sort "pks" in descending order; mark the indexies 
    [pks_rvc,locs_rvc] = findpeaks(pxx_rvc, f_rvc);            
    [pks_maxs_rvc, idx_maxs_rvc] = sort(pks_rvc, 'descend'); 

    [pxx_v,f_v] = plomb(VEL(:,2), MJD - min(MJD), 0.1, 100, 'normalized');
    [pmax_v,lmax_v] = max(pxx_v);
    f0_v = f(lmax_v);
    
    h = figure;
        hold on
%         if find(order==order_even)
%             ylim([-10 10])
%         else
%             ylim([-14 22])
%         end
        plot(f_v, pxx_v, 'k--')
        plot(f, pxx, 'r')
        plot(f_rvc, -pxx_rvc, 'b')
        legend('RV', 'Rest frame', 'Observed frame', 'Location', 'Best')
        
        for i = 1:5
            x = locs(idx_maxs(i));  % locations of the largest peaks -> harmonics
            y = pks_maxs(i);
            text(x, y, ['\leftarrow', num2str(1/x, '%3.2f')]);
            
            x_rvc = locs_rvc(idx_maxs_rvc(i));  % locations of the largest peaks -> harmonics
            y_rvc = pks_maxs_rvc(i);
            text(x_rvc, -y_rvc, ['\leftarrow', num2str(1/x_rvc, '%3.2f')]);            
        end
        
        xlabel('Frequency')
        ylabel('Normalized Power')
        title_name = ['Order', num2str(order)];
        title(title_name);
        hold off

    out_eps = [title_name, '.eps'];
    print(out_eps, '-depsc')
    close(h);
end