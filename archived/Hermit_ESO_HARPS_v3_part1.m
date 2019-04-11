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
star        = 'Gl479';
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