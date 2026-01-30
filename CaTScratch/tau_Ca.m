function tau = tau_Ca(x, dias, amp, t, Ca_paced, startpace, sys_locs, endpace, ttp, ii, i, maxcells_per_plot, p)
ca_val = dias + x.*amp;
tau_index1 = find(Ca_paced(startpace+sys_locs-1:endpace)<=ca_val,1,'first');

if isempty(tau_index1) == 0
    tau_index2 = tau_index1 - 1;
    tau_i1 = t(tau_index1);
    tau_i2 = t(tau_index2);
    ca_i1 = Ca_paced(tau_index1+startpace+sys_locs-2);
    ca_i2 = Ca_paced(tau_index2+startpace+sys_locs-2);
    ca_interp = [ca_i1; ca_i2];
    t_interp = [tau_i1; tau_i2];
    tau = interp1(ca_interp,t_interp,ca_val)+ttp;
    JM_subplot_inloop(tau+t(startpace), interp1(t,Ca_paced,tau+t(startpace)),'Time (s)','F/F0',fignum(3,i,maxcells_per_plot),p,i,2,maxcells_per_plot,'Transient Parameters');
    hold on;
else
    disp(['Peak ',num2str(ii), ' of cell ',num2str(i), ' unable to find tau.'])
    tau = NaN;
end
end
