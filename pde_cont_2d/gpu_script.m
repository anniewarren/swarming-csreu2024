%           no  n     eps     
runs = [    1   128   0.1;
            2   256   0.01;
            3   512   0.001;
            4   1024   0.0001;
            %5   1024  0.00625;
            %6   1024  0.001;
            %7   2048  0.0005;
            %8   2048  0.0001;
        ];

for i = 1:length(runs)
    name = [num2str(i,'%02u') 'n' num2str(runs(i,2),'%04u') 'eps' num2str(1000*runs(i,3),'%04u')];
    mu_list_name = genvarname([name 'mulist']);
    vac_list_name = genvarname([name 'vaclist']);
    [mu_list_out, vac_list_out] = gpu_pde_cont_2d(runs(i,2),runs(i,3),true);
    disp(['completed run ' num2str(i) ', saving data...'])
    eval([mu_list_name '= mu_list_out;']);
    eval([vac_list_name '= vac_list_out;']);
    if i == 1
        save("/export/scratch/amwarscratch/gpuRunsData",mu_list_name,vac_list_name);
    else
        save("/export/scratch/amwarscratch/gpuRunsData",mu_list_name,vac_list_name,"-append");
    end    
end    
