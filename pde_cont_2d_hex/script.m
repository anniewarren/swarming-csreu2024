%           no  n     eps     
runs = [    1   128   0.1;
            2   128   0.01;
            3   256   0.005;
            %4   512   0.0025;
        ];

for i = 1:length(runs)
    name = [num2str(runs(i,1),'%02u') 'n' num2str(runs(i,2),'%04u') 'eps' num2str(100000*runs(i,3),'%05u')];
    mu_list_name = genvarname([name 'mulist']);
    vac_list_name = genvarname([name 'vaclist']);
    [mu_list_out, vac_list_out] = pde_cont_2d_hex_run(runs(i,2),runs(i,3),false,false,false);
    disp(['completed run ' num2str(i) ', saving data...'])
    eval([mu_list_name '= mu_list_out;']);
    eval([vac_list_name '= vac_list_out;']);
    if runs(i,1) == 1
        save("C:\Users\amwar\swarming\figures\hex_vac_scaling_try1.mat",mu_list_name,vac_list_name);
        %save("/export/scratch/amwarscratch/hex_vac_scaling",mu_list_name,vac_list_name);
    else
        save("C:\Users\amwar\swarming\figures\hex_vac_scaling_try1.mat",mu_list_name,vac_list_name,"-append");
        %save("/export/scratch/amwarscratch/hex_vac_scaling",mu_list_name,vac_list_name,"-append");
    end    
end    