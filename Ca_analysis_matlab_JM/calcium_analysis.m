function calcium_analysis(tablename,savename,fps,n_paces,pacingfreq,folder,pk_threshold,paceability_threshold)

% define parameters
maxcells_per_plot=12;
remove_prepace_percent = 0.1; %intentional offset
prepace_time = round(remove_prepace_percent*n_paces+1);
Ca_table = readtable(tablename);
Ca_init = table2cell(Ca_table);
loop_length = size(Ca_init,2);
paceable_index = 0;
paceable_cell = [];
unpaceable_index = 0;
unpaceable_cell = [];

p = numSubplots(min([maxcells_per_plot, loop_length]));

for i = 1:loop_length

    try
        clearvars Ca_inter Ca_initial time_i Ca time

        Ca_inter = cell2mat(Ca_init(2:end,i));
        Ca_initial = Ca_inter(~isnan(Ca_inter));
        time_i = (0:1:length(Ca_initial)-1)./fps;
        time = time_i;
        Ca = Ca_initial;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Normalization + Filter
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % Plot raw + filtered data
        JM_subplot_inloop(time,Ca,'Time (s)','Raw F',fignum(1,i,maxcells_per_plot),p,i,0,maxcells_per_plot,'Raw + Filtered F');
        %Ca = sgolayfilt(wdenoise(medfilt1(Ca,3)),4,5);
        Ca = smooth(Ca, 15, 'rloess')
        JM_subplot_inloop(time,Ca,'Time (s)','Raw F',fignum(1,i,maxcells_per_plot),p,i,0,maxcells_per_plot,'Raw + Filtered F');

        % Identify baseline points
        F_zero = movmin(Ca,round(fps/pacingfreq*2));

        % Normalization
        norm_Ca = (Ca-F_zero)./F_zero;
        diffnorm_Ca_i = [0; diff(norm_Ca)];
        [~,diffpk_loc_i] = findpeaks(diffnorm_Ca_i,'MinPeakProminence',max(diffnorm_Ca_i)/4,'MinPeakDistance',fps/pacingfreq*(1-paceability_threshold));
        shift = diffpk_loc_i(prepace_time)-round(0.05*fps);
        norm_Ca = norm_Ca(shift:end);
        time = time(shift:end)-time(shift);
        diffnorm_Ca = [0; diff(norm_Ca)];
        [dv,dv_locs] = findpeaks(-diffnorm_Ca,'MinPeakProminence',max(-diffnorm_Ca)/4,'MinPeakDistance',fps/pacingfreq*(1-paceability_threshold));
        dv = dv.*fps;
        [uv,uv_locs] = findpeaks(diffnorm_Ca,'MinPeakProminence',max(diffnorm_Ca)/4,'MinPeakDistance',fps/pacingfreq*(1-paceability_threshold));
        uv = uv.*fps;

        % Plot normalized data
        JM_subplot_inloop(time,norm_Ca,'Time (s)','F/F0',fignum(2,i,maxcells_per_plot),p,i,0,maxcells_per_plot,'Normalized After Removing Photobleaching Effect');

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Check for paceable cells
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        beating_freqs = check_paceable(time,n_paces,pacingfreq,norm_Ca,fps);

        if beating_freqs >= (1-paceability_threshold)*pacingfreq && beating_freqs <= (1+paceability_threshold)*pacingfreq
            binary_paceable(i) = 1;
            paceable_index = paceable_index + 1;
            paceable_cell(paceable_index) = i;
            figure(fignum(2,i,maxcells_per_plot));
            subplot(p(1),p(2),1+mod(i-1,maxcells_per_plot));
            title('PACEABLE','Color','green')
        else
            binary_paceable(i) = 0;
            unpaceable_index = unpaceable_index + 1;
            unpaceable_cell(unpaceable_index) = i;
            figure(fignum(2,i,maxcells_per_plot));
            subplot(p(1),p(2),1+mod(i-1,maxcells_per_plot));
            title('UNPACEABLE','Color','red')
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Quantify Calcium Transient Parameters
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % Focus only on paced portion of data
        pacingtime = (n_paces-prepace_time+1)/pacingfreq;
        t = time(time<=pacingtime);
        Ca_paced = norm_Ca(time<=pacingtime);
        paced_uv = uv(time(uv_locs)<=pacingtime);
        paced_uv_locs = uv_locs(time(uv_locs)<=pacingtime);
        paced_dv = dv(time(dv_locs)<=pacingtime);
        paced_dv_locs = dv_locs(time(dv_locs)<=pacingtime);
        JM_subplot_inloop(t, Ca_paced,'Time (s)','F/F0',fignum(3,i,maxcells_per_plot),p,i,0,maxcells_per_plot,'Transient Parameters');
        hold on;

        if ismember(i,paceable_cell) == 1
            for ii = 1:length(paced_uv)-1
                startpace = paced_uv_locs(ii);
                endpace = paced_uv_locs(ii+1);
                [sys(ii), sys_locs(ii)] = max(Ca_paced(startpace:endpace));
                sys_locs_ppe(ii) = startpace+sys_locs(ii)-1;
                JM_subplot_inloop(t(startpace+sys_locs(ii)-1), sys(ii),'Time (s)','F/F0',fignum(3,i,maxcells_per_plot),p,i,2,maxcells_per_plot,'Transient Parameters');
                hold on;
                [dias(ii), dias_locs(ii)] = min(Ca_paced(startpace:endpace));
                JM_subplot_inloop(t(startpace+dias_locs(ii)-1), dias(ii),'Time (s)','F/F0',fignum(3,i,maxcells_per_plot),p,i,2,maxcells_per_plot,'Transient Parameters');
                hold on;
                amp(ii) = sys(ii)-dias(ii);
                ttp(ii) = t(sys_locs(ii));
                tau50(ii) = tau_Ca(0.5, dias(ii), amp(ii), t, Ca_paced, startpace, sys_locs(ii), endpace, ttp(ii), ii, i, maxcells_per_plot, p);
                tau90(ii) = tau_Ca(0.1, dias(ii), amp(ii), t, Ca_paced, startpace, sys_locs(ii), endpace, ttp(ii), ii, i, maxcells_per_plot, p);
                x = (t(paced_uv_locs(ii)-1)):1/fps:t(paced_uv_locs(ii)+1);
                y = uv(ii).*(x-t(paced_uv_locs(ii)))+Ca_paced(paced_uv_locs(ii));
                JM_subplot_inloop(x,y,'Time (s)','F/F0',fignum(3,i,maxcells_per_plot),p,i,0,maxcells_per_plot,'Transient Parameters');
                hold on;
            end
            for ii=1:length(paced_dv)-1
                x = (t(paced_dv_locs(ii)-1)):1/fps:t(paced_dv_locs(ii)+1);
                y = -dv(ii).*(x-t(paced_dv_locs(ii)))+Ca_paced(paced_dv_locs(ii));
                JM_subplot_inloop(x,y,'Time (s)','F/F0',fignum(3,i,maxcells_per_plot),p,i,0,maxcells_per_plot,'Transient Parameters');
                hold on;
            end

            Amp_mean(i) = mean(amp); % mean amplitude
            Amp_std(i) = std(amp); % std amplitude
            Dias_mean(i) = mean(dias); % mean amplitude
            Dias_std(i) = std(dias); % std amplitude

            %%% upstroke velocity
            UV_mean(i) = mean(uv); % mean upstroke velocity
            UV_std(i) = std(uv); % std upstroke velocity

            T50_mean(i) = mean(tau50(tau50<(1/pacingfreq))); % mean T50
            T90_mean(i) = mean(tau90(tau90<(1/pacingfreq))); % mean T90
            T50_std(i) = std(tau50(tau50<(1/pacingfreq))); % std T50
            T90_std(i) = std(tau90(tau90<(1/pacingfreq))); % std T90
            ttp_mean(i) = mean(ttp); % mean T90
            ttp_std(i) = std(ttp); % std T50
            DV_mean(i) = mean(dv); % mean downstroke velocity
            DV_std(i) = std(dv); % std downstroke velocity

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Quantify Number/Duration of aCREs
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            [sys(end+1), sys_locs(end+1)] = max(Ca_paced(paced_uv_locs(end):end));
            sys_locs_ppe(end+1) = sys_locs(end) + paced_uv_locs(end);
            [acre_potential,acre_potential_loc] = findpeaks(norm_Ca,'MinPeakProminence',Amp_mean(i).*pk_threshold);
            JM_subplot_inloop(time, norm_Ca,'Time (s)','F/F0',fignum(4,i,maxcells_per_plot),p,i,0,maxcells_per_plot,'aCREs');
            hold on;
            JM_subplot_inloop(time(acre_potential_loc),acre_potential,'Time (s)','F/F0',fignum(4,i,maxcells_per_plot),p,i,1,maxcells_per_plot,'aCREs');
            hold on;
            acre_idx_pre = 0;
            acre_idx_post = 0;
            acre_pre = 0;
            acre_post = 0;
            for ii = 1:length(acre_potential_loc)
                if any(abs(sys_locs_ppe'-acre_potential_loc(ii))<=round(0.05*fps))
                    JM_subplot_inloop_text(time(acre_potential_loc(ii)),acre_potential(ii),'N',fignum(4,i,maxcells_per_plot),p,i,maxcells_per_plot,'aCREs');
                    hold on;
                else
                    if acre_potential_loc(ii)>(sys_locs_ppe(end)+round(0.025*fps))
                        acre_idx_post = acre_idx_post+1;
                        acre_post(acre_idx_post) = acre_potential_loc(ii);
                        JM_subplot_inloop_text(time(acre_potential_loc(ii)),acre_potential(ii),'A',fignum(4,i,maxcells_per_plot),p,i,maxcells_per_plot,'aCREs');
                    else
                        acre_idx_pre = acre_idx_pre+1;
                        acre_pre(acre_idx_pre) = acre_potential_loc(ii);
                        JM_subplot_inloop_text(time(acre_potential_loc(ii)),acre_potential(ii),'A',fignum(4,i,maxcells_per_plot),p,i,maxcells_per_plot,'aCREs');
                    end
                end
            end
            acre_pacing(i) = length(acre_pre)-1;
            acre_ppe(i) = length(acre_post)-1;
            acre_tot(i) = acre_pacing(i)+acre_ppe(i);
        else
            JM_subplot_inloop(NaN, NaN,'Time (s)','Unpaceable',fignum(3,i,maxcells_per_plot),p,i,0,maxcells_per_plot,'Transient Parameters');
            Amp_mean(i) = NaN;
            Amp_std(i) = NaN;
            Dias_mean(i) = NaN;
            Dias_std(i) = NaN;
            UV_mean(i) = NaN;
            UV_std(i) = NaN;
            DV_mean(i) = NaN;
            DV_std(i) = NaN;
            T50_mean(i) = NaN;
            T90_mean(i) = NaN;
            T50_std(i) = NaN;
            T90_std(i) = NaN;
            ttp_mean(i) = NaN;
            ttp_std(i) = NaN;
            acre_pacing(i) = NaN;
            acre_ppe(i) = NaN;
            acre_tot(i) = NaN;
        end

    catch
        warning(['Analysis of Cell',num2str(i),' skipped due to error']);

        JM_subplot_inloop(NaN, NaN,'Time (s)','Unpaceable',fignum(3,i,maxcells_per_plot),p,i,0,maxcells_per_plot,'Transient Parameters');
        JM_subplot_inloop(NaN, NaN,'Time (s)','Unpaceable',fignum(4,i,maxcells_per_plot),p,i,0,maxcells_per_plot,'aCREs');

        Amp_mean(i) = NaN;
        Amp_std(i) = NaN;
        Dias_mean(i) = NaN;
        Dias_std(i) = NaN;
        UV_mean(i) = NaN;
        UV_std(i) = NaN;
        DV_mean(i) = NaN;
        DV_std(i) = NaN;
        T50_mean(i) = NaN;
        T90_mean(i) = NaN;
        T50_std(i) = NaN;
        T90_std(i) = NaN;
        ttp_mean(i) = NaN;
        ttp_std(i) = NaN;
        acre_pacing(i) = NaN;
        acre_ppe(i) = NaN;
        acre_tot(i) = NaN;
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 6: Output to CSV File
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

varNames = {'Cell_Number','Amp_mean_F_F0','Amp_std_F_F0','Dias_mean_F_F0','Dias_std_F_F0','UV_mean_F_F0_s','UV_std_F_F0_s',...
    'DV_mean_F_F0_s','DV_std_F_F0_s','T50_mean_s','T50_std_s','T90_mean_s','T90_std_s','TTP_mean_s','TTP_std_s','aCRE_pacing','aCRE_PPE','aCRE_total'};

relabel_table = readtable(tablename,'ReadVariableNames',0);

for names = 1:loop_length
    cell_name(names) = relabel_table{1,names};
end

T = table(cell_name',Amp_mean',Amp_std',Dias_mean',Dias_std',...
    UV_mean',UV_std',DV_mean',DV_std',...
    T50_mean',T50_std',T90_mean',T90_std',ttp_mean', ttp_std', acre_pacing',acre_ppe',acre_tot',...
    'VariableNames',varNames);

T2 = rmmissing(T);

writetable(T2,[folder, '/', savename,'_Results.csv'],'WriteVariableNames',1)

varNames_real = {'Cell_Name','Mean Amplitude (F/F0)','Std Amplitude (F/F0)','Mean Diastolic (F/F0)','Std Diastolic (F/F0)',...
    'Mean Upstroke Velocity (F/F0/s)','Std Upstroke Velocity (F/F0/s)',...
    'Mean Downstroke Velocity (F/F0/s)','Std Downstroke Velocity (F/F0/s)',...
    'Mean T50 (s)','Std T50 (s)','Mean T90 (s)','Std T90 (s)','Mean TTP (s)','Std TTP (s)','aCRE_pacing (n)','aCRE_PPE (n)','aCRE_total (n)'};

FigureFileName = strcat(folder, '/', savename, '_allfigures.fig');
FigList = findobj(allchild(0), 'flat', 'Type', 'figure'); savefig(FigList,FigureFileName);

end
