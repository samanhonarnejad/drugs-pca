% Read, Organize and save data for the experiment 
% Population mean data for each cell line
% Single cell data for each 

clear;

cell_lines = {
    '10a'
    };
timepoints = {
    '24h'
    '48h'
    };
drugs = {
    'AZD6244'
    'BEZ235'
    'Erotinib'
    'Gefitinib'
    'Lapatinib'
    'MK2206'
    'PD0325901'
    'PP242'
    'Triciribine'
    };
signals_647 = { % 18 signals from Alexa647 channels
    'pRSK_Ser380'
    'pERK_Thr202_Tyr204'
    'pAKT_Ser473'
    'Foxo3a'
    'pS6_Ser235_Ser236'
    'p4EBP1_Thr37_46'
    'pCDK2_Tyr15'
    'pCDK1_Tyr15'
    'pP57_Thr310'
    'pP27_Ser10'
    'pP27_Thr187'
    'Survivin'
    'p27'
    'p21'
    'p57'
    'FoxM1'
    'CyclinB'
    'CyclinA'
    };

doses = [(10 * (3.1623.^(0:-1:-6)))*10^-6 0]; % drug concentration in M
doses = doses(8:-1:1);


% Handling single-cell data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear cell_counter;
for cell_line_id = 1:1
    for drug_id = 1:9
        for timepoint_id = 1:2
            for plate_id = 1:3
                for i = 'A':'H'
                    for j = 1:12
                        foldername = sprintf('%s_%s_%s_P%d/',char(cell_lines(cell_line_id)),char(drugs(drug_id)),char(timepoints(timepoint_id)),plate_id);
                        L = dir([foldername '*result.' i num2str(j) '[' '*Selected[1].csv']);
                        if length(L) > 0
                            row = i - 64;
                            col = j;
                            
                            signals_647_id = (plate_id-1)*6 + floor((col-1)/2)+1;
                            dose_id = row;
                            rep = rem(col-1,2)+1;

                            
                            filename = sprintf('%s_%s_%s_P%d/%s',char(cell_lines(cell_line_id)),char(drugs(drug_id)),char(timepoints(timepoint_id)),plate_id,L.name);
                            csv_content = dlmread(filename,',',1,16);
                            [n p] = size(csv_content);
                            
                            cell_counter(cell_line_id,drug_id,timepoint_id,dose_id,signals_647_id,rep) = n; % cell number for each well to be used later
                            
                            for cell_id = 1:n
                                
                                single_cell_DNA_content_data(cell_line_id,drug_id,timepoint_id,dose_id,signals_647_id,rep,cell_id) = csv_content(cell_id,6);
                                
                                single_cell_Alexa647_nuc_cyt_ratio_data(cell_line_id,drug_id,timepoint_id,dose_id,signals_647_id,rep,cell_id) = csv_content(cell_id,1);
                                single_cell_Alexa647_nuc_data(cell_line_id,drug_id,timepoint_id,dose_id,signals_647_id,rep,cell_id) = csv_content(cell_id,3);
                                single_cell_Alexa647_cell_data(cell_line_id,drug_id,timepoint_id,dose_id,signals_647_id,rep,cell_id) = csv_content(cell_id,4);
                                
                                single_cell_Alexa568_nuc_data(cell_line_id,drug_id,timepoint_id,dose_id,signals_647_id,rep,cell_id) = csv_content(cell_id,5);
                                
                                single_cell_Alexa568_nuc_data(cell_line_id,drug_id,timepoint_id,dose_id,signals_647_id,rep,cell_id) = csv_content(cell_id,5);
                            end
                        end
                    end
                end
            end
        end
    end
end

max_cell_number = max(reshape(cell_counter,[],1));
for cell_line_id = 1:1
    for drug_id = 1:9
        for timepoint_id = 1:2
            for signals_647_id = 1:18
                for dose_id = 1:8
                    for rep = 1:2
                        if cell_counter(cell_line_id,drug_id,timepoint_id,dose_id,signals_647_id,rep) < max_cell_number
                            for cell_id = cell_counter(cell_line_id,drug_id,timepoint_id,dose_id,signals_647_id,rep)+1:max_cell_number
                                
                                single_cell_DNA_content_data(cell_line_id,drug_id,timepoint_id,dose_id,signals_647_id,rep,cell_id) = nan;
                                
                                single_cell_Alexa647_nuc_cyt_ratio_data(cell_line_id,drug_id,timepoint_id,dose_id,signals_647_id,rep,cell_id) = nan;
                                single_cell_Alexa647_nuc_data(cell_line_id,drug_id,timepoint_id,dose_id,signals_647_id,rep,cell_id) = nan;
                                single_cell_Alexa647_cell_data(cell_line_id,drug_id,timepoint_id,dose_id,signals_647_id,rep,cell_id) = nan;
                                
                                single_cell_Alexa568_nuc_data(cell_line_id,drug_id,timepoint_id,dose_id,signals_647_id,rep,cell_id) = nan;
                            end
                        end
                    end
                end
            end
        end
    end
end

% Handling population average data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for cell_line_id = 1:1
    for drug_id = 1:9
        for timepoint_id = 1:2
            for signals_647_id = 1:18
                for dose_id = 1:8
                    for rep = 1:2
                        Alexa647_nuc_cyt_ratio_data(cell_line_id,drug_id,timepoint_id,dose_id,signals_647_id,rep) = nanmean(squeeze(single_cell_Alexa647_nuc_cyt_ratio_data(cell_line_id,drug_id,timepoint_id,dose_id,signals_647_id,rep,:)));
                        Alexa647_nuc_data(cell_line_id,drug_id,timepoint_id,dose_id,signals_647_id,rep) = nanmean(squeeze(single_cell_Alexa647_nuc_data(cell_line_id,drug_id,timepoint_id,dose_id,signals_647_id,rep,:)));
                        Alexa647_cell_data(cell_line_id,drug_id,timepoint_id,dose_id,signals_647_id,rep) = nanmean(squeeze(single_cell_Alexa647_cell_data(cell_line_id,drug_id,timepoint_id,dose_id,signals_647_id,rep,:)));
                        
                        Alexa568_nuc_data(cell_line_id,drug_id,timepoint_id,dose_id,signals_647_id,rep) = nanmean(squeeze(single_cell_Alexa568_nuc_data(cell_line_id,drug_id,timepoint_id,dose_id,signals_647_id,rep,:)));
                    end
                end
            end
        end
    end
end

% Normalizing data to control (untreated) condition
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for cell_line_id = 1:1
    for drug_id = 1:9
        for timepoint_id = 1:2
            for signals_647_id = 1:18
                for dose_id = 1:8
                    for rep = 1:2
                        normalized_Alexa647_nuc_cyt_ratio_data(cell_line_id,drug_id,timepoint_id,dose_id,signals_647_id,rep) = ...
                            Alexa647_nuc_cyt_ratio_data(cell_line_id,drug_id,timepoint_id,dose_id,signals_647_id,rep) / mean(Alexa647_nuc_cyt_ratio_data(cell_line_id,drug_id,timepoint_id,1,signals_647_id,:));
                        normalized_Alexa647_nuc_data(cell_line_id,drug_id,timepoint_id,dose_id,signals_647_id,rep) = ...
                            Alexa647_nuc_data(cell_line_id,drug_id,timepoint_id,dose_id,signals_647_id,rep) / mean(Alexa647_nuc_data(cell_line_id,drug_id,timepoint_id,1,signals_647_id,:));
                        normalized_Alexa647_cell_data(cell_line_id,drug_id,timepoint_id,dose_id,signals_647_id,rep) = ...
                            Alexa647_cell_data(cell_line_id,drug_id,timepoint_id,dose_id,signals_647_id,rep) / mean(Alexa647_cell_data(cell_line_id,drug_id,timepoint_id,1,signals_647_id,:));
                        normalized_Alexa568_nuc_data(cell_line_id,drug_id,timepoint_id,dose_id,signals_647_id,rep) = ...
                            Alexa568_nuc_data(cell_line_id,drug_id,timepoint_id,dose_id,signals_647_id,rep) / mean(Alexa568_nuc_data(cell_line_id,drug_id,timepoint_id,1,signals_647_id,:));
                    end
                end
            end
        end
    end
end



          
% Save data onto a ".mat" file for later use
% Note that if you're using home server the file might not save there, in
% that case use a local location like Desktop for saving the file and then
% move it

save Organize_10a_Drug_Response_Single_Cell_Data.mat;


