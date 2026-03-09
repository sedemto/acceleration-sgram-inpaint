function [solution] = inpaint(sgram,param,paramsolver)

mask = sgram(1,:) ~= 0;
insig = param.F_adj(sgram);

soft = @(z, lambda) sign(z).*max(abs(z) - lambda, 0);
param.proj = @(x) x.*(1-mask) + sgram.*mask;
param.prox = @(z) soft(z, paramsolver.lambda/paramsolver.sigma);


    
paramsolver.x0 = insig;
paramsolver.y0 = zeros(size(param.F(insig)));
paramsolver.z0 = zeros(size(param.F(insig)));

%% initialization
x = paramsolver.x0;
y = paramsolver.y0;
z = paramsolver.z0;

tau = paramsolver.tau;
sigma = paramsolver.sigma;
alpha = paramsolver.alpha;
eta = paramsolver.eta;

%% iteration of generalized Chambolle-Pock algorithm

for i = 1:paramsolver.I

    tmp_r = y + eta*param.F(x - tau*(param.F_adj(z)+param.F_adj(y)));
    r = tmp_r - eta*param.proj(tmp_r/eta);

    p = x - tau*(param.F_adj(z) + param.F_adj(r));

    tmp_q = z + sigma*param.F(2*p - x);
    q = tmp_q - sigma*param.prox(tmp_q/sigma);
    
    x = x + alpha*(p - x);
    y = y + alpha*(r - y);
    z = z + alpha*(q - z);
    
end

solution = param.proj(param.F(x));

end

