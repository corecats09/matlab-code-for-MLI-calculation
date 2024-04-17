param_log_adj = zeros(15,15);
for m = 1:15
cvx_begin
A = []; 
variable A(m, 1); 
minimize( norm(tbl_ele_adj (:,1:m) *A-mean_log, 'fro') );
subject to
A >= 0; % 所有元素非负
cvx_end; 
param_log_adj(1:m,m)=A;
end
param_sta_adj=zeros(15,15);
for m = 1:15
cvx_begin
A = []; 
variable A(m, 1); 
minimize( norm(tbl_ele_adj (:,1:m) *A-mean_sta, 'fro') );
subject to
A >= 0; % 所有元素非负
cvx_end;
param_sta_adj(1:m,m)=A;
end



for n = 1:15
residual_sta(:,n) = tbl_ele_adj*param_sta_adj(:,n)-mean_sta;
std_residual_sta(1,n) = std (residual_sta(:,n));
end
for n = 1:15
lnLmax_sta(:,n) = length(mean_sta).*log(sqrt(2.*pi).*std_residual_sta(1,n))- norm(residual_sta(:,n)-mean(residual_sta(:,n)))^4./2./(std_residual_sta(1,n)^2);
end
for n = 1:15
cov_temp =std_residual_sta(1,n)^2* inv(tbl_ele_adj(:,1:n)'* tbl_ele_adj(:,1:n));
det_cov_sta(1,n) = det(cov_temp);
end
for n = 1:15
lnMLI_sta(1,n) = n/2.*log(2*pi)+lnLmax_sta(1,n)+0.5.*log(det_cov_sta(1,n))-n.*log(1);
end

for n = 1:15
residual_log(:,n) = tbl_ele_adj*param_log_adj(:,n)-mean_log;
std_residual_log(1,n) = std (residual_log(:,n));
end
for n = 1:15
lnLmax_log(:,n) = length(mean_log).*log(sqrt(2.*pi).*std_residual_log(1,n))-norm(residual_log(:,n)-mean(residual_log(:,n)))^4./2./(std_residual_log(1,n)^2);
end
for n = 1:15
cov_temp =std_residual_log(1,n)^2  * inv(tbl_ele_adj(:,1:n)'* tbl_ele_adj(:,1:n));
det_cov_log(1,n) = det(cov_temp);
end
for n = 1:15
lnMLI_log(1,n) = n/2.*log(2*pi)+lnLmax_log(1,n)+0.5.*log(det_cov_log(1,n))-n.*log(1);
end

scatter (linspace(1,15,15),lnMLI_log)
hold on
scatter (linspace(1,15,15),lnMLI_sta)
