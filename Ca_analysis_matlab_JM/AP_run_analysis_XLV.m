function AP_run_analysis(mainfile,parameterfile,folder,pk_threshold,paceability)
fprintf(parameterfile);
parameters = csvread(parameterfile); % call parameters
[~, savename, ~] = fileparts(mainfile); % naming of output excel file
AP_analysis_XLV(mainfile,savename,parameters(1,1),parameters(1,2),parameters(1,3),folder,pk_threshold,paceability);

end
