function [n_ideal] = Get_Estimated_n(cstarcp,modeledcp,nall)


cp_all = modeledcp; % so
cstarcp = cstarcp(:);
%% calculate square difference between cp_all and cstarcp

cstar = repmat(cstarcp,1,length(cp_all(1,:)));
sq = (cstar-cp_all).^2; sq=sq';
[minimized_cp ind]= nanmin(sq); n_ideal=nall(ind);n_ideal(n_ideal==1.001)=NaN;

%hist(minimized_cp)% ok, kjust checking that values are very low. about 10^-9.
ind = find(minimized_cp>0.0001); n_ideal(ind)=NaN; minimized_cp(ind)=NaN;

n_ideal = n_ideal(:);

figure
yyaxis left
plot(cstar,'r');
ylabel('CSTAR c_p [m^-^1]')
yyaxis right
plot(n_ideal,'b');
ylabel('relative real refractive index')
legend('c_p','n+0.0003i')
ax=gca;
ax.YAxis(1).Color='r';
ax.YAxis(2).Color='b';
datetick
xlabel('Time (2017)')

figure
plot(cstarcp,n_ideal,'.')

end


