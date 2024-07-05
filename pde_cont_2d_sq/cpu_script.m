%           no  n     eps     
runs = [    %1   128   0.1;
            2   128   0.01;
            3   256   0.001;
            4   512   0.0001;
        ];

for i = 1:length(runs)
    name = [num2str(i,'%02u') 'n' num2str(runs(i,2),'%04u') 'eps' num2str(100000*runs(i,3),'%05u')];
    mu_list_name = genvarname([name 'mulist']);
    vac_list_name = genvarname([name 'vaclist']);
    inf_list_name = genvarname([name 'inflist']);
    [mu_list_out, vac_list_out, inf_list_out] = cpu_pde_cont_2d(runs(i,2),runs(i,3));
    disp(['completed run ' num2str(i) ', saving data...'])
    eval([mu_list_name '= mu_list_out;']);
    eval([vac_list_name '= vac_list_out;']);
    eval([inf_list_name '= inf_list_out;']);
    if runs(i,1) == 2
        save("C:\Users\amwar\swarming\data\bif_diag_inset_cont.mat",mu_list_name,vac_list_name,inf_list_name);
        %save("/export/scratch/amwarscratch/cpuRunsData",mu_list_name,vac_list_name);
    else
        %save("C:\Users\amwar\swarming\data\bif_diag_inset_cont.mat",mu_list_name,vac_list_name,inf_list_name,"-append");
        %save("/export/scratch/amwarscratch/cpuRunsData",mu_list_name,vac_list_name,"-append");
    end    
end    