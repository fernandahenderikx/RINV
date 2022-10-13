% modified to include new Kij matrices, Jan 2021
function get_pscat(all_min_time,all_min_data)

% change to folder containing cruise of interest
cd ../
hdir = cd;

mdirL   = ([hdir '\mfiles\']); 
figpath = [hdir '\figures\'];
datapath = [hdir '\mat_files\'];
calpath  = [hdir '\instrument_files\'];

ringareafile = [calpath 'ringarea_1421.asc']; %ring area same for 1285 and 1421
dias         = load([calpath 'Dias32_b.asc']);
fact_zsc     = load([calpath 'factory_zsc_1421_July2016.asc']);

%% variables
% LISST S/N 1421, 2016 Cal prior to sikuliak
%         VCC = 33720;    % from InstrumentData.txt, changes with factory calibration
%         flref = 873.45; % from element 36 in factory_zsc file

%% process data

 % matlab date 
      sdate = datenum('1/0/2017') + all_min_time; % time in GMT hopefully
      dvec = datevec(sdate);
          
%  plot dvec and data to check timing of fsw; correct for outliers
    figure         
    plot(dvec(:,5),all_min_data(:,20),'.')
    
% there are times with all filter and times with all total           
%  fsw comes in 0-10; time looks right with some outliers
%  filter oout bad fsw sections
    cd(mdirL)
  % perform inversion for iterative sections of data using 1-10min fsw  
        vd(1:size(all_min_data,1),1:32)= NaN;
        tau_all(1:size(all_min_data,1),1:32)= NaN;
      
        zscfile=[cd '\zsc.asc'];
        logfile=[cd '\rdata.log'];  
        last_fsw = [17:18];               % first decent fsw section, First on KM1709 17:18
        %% how to select a fsw segment, enter this line and manually pick a fsw section
%         plot(all_min_data(:,25))
%               
        for mo= nanmin(dvec(:,2)):nanmax(dvec(:,2)) %month. 
            mm = find(dvec(:,2)==mo);
        for d=nanmin(dvec(mm,3)):nanmax(dvec(mm,3)) %days
            for h=0:23
              whole_seg = find(dvec(:,2)==mo & dvec(:,3)==d & dvec(:,4)==h &  dvec(:,5)>11);
              fsw_seg = find(dvec(:,2)==mo & dvec(:,3)==d & dvec(:,4)==h &  dvec(:,5)<8 & dvec(:,5)>2);
              %test to ensure whole > fsw
              ftest = nanmean(nanmean(all_min_data(whole_seg,1:32)))-nanmean(nanmean(all_min_data(fsw_seg,1:32)));
               if (isempty(fsw_seg) || size(fsw_seg,1)==1 || ftest<1)
                       fsw_seg = last_fsw;
               end
               if (isempty(whole_seg)==0)
                    last_fsw = fsw_seg; 
                    fsw=nanmean(all_min_data(fsw_seg,:));
                    %fsw=fsw.*0; % make it all zeros to test whether that improves pscat
                    part=all_min_data(whole_seg,:);
                    plot(dias,fsw(:,1:32),'k'); hold on; plot(dias,part(:,1:32),'r');% pause(1); 
                    clf
                 save([cd,'\zsc.asc'],'fsw','-ascii');
            save([cd,'\rdata.log'],'part','-ascii');
          % [scat,tau,zsc,data_out,cscat] = getscat(logfile,zscfile,1,ringareafile); %all data corrected for deepest 2m
           

             [angr,pscat,beamc,snr,ang1,ang2,rings,cscat,tau,rawdata] = readpscat(logfile, zscfile, ringareafile);

                 cscat_all(whole_seg,:) = cscat;
                 pscat_all(whole_seg,:) = pscat;


                end
            end
        end
        end
      
        
     cd(datapath)
           save('LISST_pscat.mat','sdate','vd','pscat_all','cscat_all')
  
  end
  