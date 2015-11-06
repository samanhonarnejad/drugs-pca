
clear;
load Organize_10a_Drug_Response_Single_Cell_Data.mat;

% Make excel sheets from population average data

cell_line_id = 1;
for drug_id = 1:9
    for dose_id = 1:8
        row = dose_id;
        col = 0;
        for timepoint_id = 1:2
            for signals_647_id = 1:18
                for rep = 1:2
                    col = col + 1;
                    excell_formatted_Alexa647_nuc_cyt_ratio_data(row,col) = Alexa647_nuc_cyt_ratio_data(cell_line_id,drug_id,timepoint_id,dose_id,signals_647_id,rep);
                    excell_formatted_Alexa647_nuc_data(row,col) = Alexa647_nuc_data(cell_line_id,drug_id,timepoint_id,dose_id,signals_647_id,rep);
                    excell_formatted_Alexa647_cell_data(row,col) = Alexa647_cell_data(cell_line_id,drug_id,timepoint_id,dose_id,signals_647_id,rep);
                    excell_formatted_Alexa568_nuc_data(row,col) = Alexa568_nuc_data(cell_line_id,drug_id,timepoint_id,dose_id,signals_647_id,rep);
                    
                    
                    normalized_excell_formatted_Alexa647_nuc_cyt_ratio_data(row,col) = normalized_Alexa647_nuc_cyt_ratio_data(cell_line_id,drug_id,timepoint_id,dose_id,signals_647_id,rep);
                    normalized_excell_formatted_Alexa647_nuc_data(row,col) = normalized_Alexa647_nuc_data(cell_line_id,drug_id,timepoint_id,dose_id,signals_647_id,rep);
                    normalized_excell_formatted_Alexa647_cell_data(row,col) = normalized_Alexa647_cell_data(cell_line_id,drug_id,timepoint_id,dose_id,signals_647_id,rep);
                    normalized_excell_formatted_Alexa568_nuc_data(row,col) = normalized_Alexa568_nuc_data(cell_line_id,drug_id,timepoint_id,dose_id,signals_647_id,rep);

                end
            end
        end
    end
    filename1 = sprintf('Mean_Analysis_10a/%s_Alexa647_nuc_cyt_ratio_data.csv',char(drugs(drug_id)));
    filename2 = sprintf('Mean_Analysis_10a/%s_Alexa647_nuc_data.csv',char(drugs(drug_id)));
    filename3 = sprintf('Mean_Analysis_10a/%s_Alexa647_cell_data.csv',char(drugs(drug_id)));
    filename4 = sprintf('Mean_Analysis_10a/%s_Alexa568_nuc_data.csv',char(drugs(drug_id)));
    dlmwrite(filename1,excell_formatted_Alexa647_nuc_cyt_ratio_data);
    dlmwrite(filename2,excell_formatted_Alexa647_nuc_data);
    dlmwrite(filename3,excell_formatted_Alexa647_cell_data);
    dlmwrite(filename4,excell_formatted_Alexa568_nuc_data);
    
    filename1 = sprintf('Mean_Analysis_10a/%s_normalized_Alexa647_nuc_cyt_ratio_data.csv',char(drugs(drug_id)));
    filename2 = sprintf('Mean_Analysis_10a/%s_normalized_Alexa647_nuc_data.csv',char(drugs(drug_id)));
    filename3 = sprintf('Mean_Analysis_10a/%s_normalized_Alexa647_cell_data.csv',char(drugs(drug_id)));
    filename4 = sprintf('Mean_Analysis_10a/%s_normalized_Alexa568_nuc_data.csv',char(drugs(drug_id)));
    dlmwrite(filename1,normalized_excell_formatted_Alexa647_nuc_cyt_ratio_data);
    dlmwrite(filename2,normalized_excell_formatted_Alexa647_nuc_data);
    dlmwrite(filename3,normalized_excell_formatted_Alexa647_cell_data);
    dlmwrite(filename4,normalized_excell_formatted_Alexa568_nuc_data);
    
    
end
                