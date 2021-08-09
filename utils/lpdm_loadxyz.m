addpath('/n/home03/xwei/toolbox');

fn='./OUT_LPDM/';     %%% change to your LPDM output directory
fn_save='./OUT_LPDM/';    %%% change to where you want to put the converted .mat data
caseid='RCE_192n_lpdm_example_test'; %%% change to your casename

%%%% preset parameters
%%% para.nx and para.ny are domain size, i.e., nx_gl and ny_gl
%%%
para.nx    =16;   para.ny    =16; %domain size, i.e., nx_gl and ny_gl
para.lpdm_num =1200; %lpdm_num in lpdm_mod.f90
para.nsub=4; %lpdm_sub in lpdm_mod.f90
seg = para.nx *para.ny *para.lpdm_num /para.nsub;
Nt= 240; ndt= 8;  % Nt is nstop in prm, ndt is lpdm_ndt in prm

for nsub= 1:para.nsub
    fnn1=[fn,num2str(caseid),'_sub',num2str(nsub),'.xyz'];
    fid1=fopen(fnn1);
    xyz=zeros(seg,Nt/ndt,'int32');

    for t=1:Nt/ndt+1
        fread(fid1,1,'int32');
        xyz(:,t) =fread(fid1,seg,'int32');
        fread(fid1,1,'int32');
    end

    fclose(fid1);
    fnn=[fn_save,num2str(caseid),'_sub',num2str(nsub),'.mat'];
    save(fnn,'xyz','-v7.3');

end
