% Branch from Hermit_coeff_nor.m

% Integrated from Hermit_coeff_NOR_0720_GI479.m

%%%%%%%%%%
% Update %
%%%%%%%%%%


%%%%%%%%%%%%%%
% Parameters %
%%%%%%%%%%%%%%
% star        = 'HD41248';
star        = 'Gl628';
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
    % b           = f.b;  % shift
    b           = RV_HARPS(n);  % shift
    % plot(v,A, v,f(v)) % test %
    RV_gauss(n) = b;    

    parfor order = 0:ORDER
        % temp                    = hermite_nor(order, v + 9.2) * grid_size;
        temp                    = hermite_nor(order, v - RVC) * grid_size;
        temp_rvc                = hermite_nor(order, v - b) * grid_size;
        coeff(order+1, n)       = sum(A .* temp);  
        coeff_rvc(order+1, n)   = sum(A .* temp_rvc); 
    end
    waitbar( n / N_FILE )
end
close(h)  

cd (['../', star, '/'])

for order = 0:ORDER
    
    [pxx,f] = plomb(coeff(order+1, :), MJD - min(MJD), 0.5,  'normalized');
    [pmax,lmax] = max(pxx);
    f0 = f(lmax);
    disp(['T_planet: ', num2str(1/f0)]);

    [pxx_rvc,f_rvc] = plomb(coeff_rvc(order+1, :), MJD - min(MJD), 0.5,  'normalized');
    [pmax_rvc,lmax_rvc] = max(pxx_rvc);
    f0_rvc = f_rvc(lmax_rvc);
    disp(['T_activity: ', num2str(1/f0_rvc)]);

    % pxx_eve = envelope(pxx,2,'peak');
    pxx_eve = pxx;
    [pks,locs] = findpeaks(pxx_eve, f);                 % find all the peaks in (pxx, f)
    [pks_maxs, idx_maxs] = sort(pks, 'descend');    % sort "pks" in descending order; mark the indexies 
    % pxx_rvc_eve = envelope(pxx_rvc,2,'peak');
    pxx_rvc_eve = pxx_rvc;
    [pks_rvc,locs_rvc] = findpeaks(pxx_rvc_eve, f_rvc);            
    [pks_maxs_rvc, idx_maxs_rvc] = sort(pks_rvc, 'descend'); 

     [pxx_v,f_v] = plomb(RV_gauss, MJD - min(MJD), 0.5, 'normalized');
     [pmax_v,lmax_v] = max(pxx_v);
     f0_v = f(lmax_v);
    
    h = figure;
%         if find(order==order_even)
%             ylim([-10 10])
%         else
%             ylim([-14 22])
%         end
        % ylim([-12 15])
        
        semilogx(1./f_v, pxx_v, 'k--', 1./f, pxx_eve, 'r', 1./f_rvc, -pxx_rvc_eve, 'b')
        xlim([2 150])
        % plot(f, pxx, 'r')
        % plot(f_rvc, -pxx_rvc, 'b')
        % semilogx(1./f, pxx_eve, 'r')
        % semilogx(1./f_rvc, -pxx_rvc_eve, 'b')
        
        legend('RV', 'Rest frame', 'Observed frame', 'Location', 'Best')
        hold on

        for i = 1:10
            x = locs(idx_maxs(i));  % locations of the largest peaks -> harmonics
            y = pks_maxs(i);
            if (1/x<150) && (1/x>2)
                text(1/x, y, ['\leftarrow', num2str(1/x, '%3.2f')], 'fontsize', 8);
            end
            
            x_rvc = locs_rvc(idx_maxs_rvc(i));  % locations of the largest peaks -> harmonics
            y_rvc = pks_maxs_rvc(i);
            if (1/x_rvc<150) && (1/x_rvc>2)
                text(1/x_rvc, -y_rvc, ['\leftarrow', num2str(1/x_rvc, '%3.2f')], 'fontsize', 8);      
            end
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