% Plot Figures from the data for experiment "Organize_10a_Drug_Response_Single_Cell_Data"

clear;
close all;
load Organize_10a_Drug_Response_Single_Cell_Data;

colors = [
    1 0 0
    0 0.50 0
    0 0 1
    0.75 0 0.75
    1 1 0
    ];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plotting scatter plot for two signals at 3 different doses
figure(1);
cell_line_id = 1;

drug_id1 = 1;
drug_id2 = 2;
drug_id3 = 3;
drug_id4 = 4;
drug_id5 = 5;
drug_id6 = 6;
drug_id7 = 7;
drug_id8 = 8;
drug_id9 = 9;

timepoint_id1 = 1;
timepoint_id2 = 2;

signals_647_id1 = 1;
signals_647_id2 = 2;
signals_647_id3 = 3;
signals_647_id4 = 4;
signals_647_id5 = 5;
signals_647_id6 = 6;
signals_647_id7 = 7;
signals_647_id8 = 8;
signals_647_id9 = 9;
signals_647_id10 = 10;
signals_647_id11 = 11;
signals_647_id12 = 12;
signals_647_id13 = 13;
signals_647_id14 = 14;
signals_647_id15 = 15;
signals_647_id16 = 16;
signals_647_id17 = 17;
signals_647_id18 = 18;


j = 0;
for dose_id = [1 5 6 8]
    j = j+1;
    x1 = reshape((single_cell_Alexa647_nuc_data(cell_line_id,drug_id8,timepoint_id1,dose_id,signals_647_id6,:,:)),[],1);
    y1 = reshape((single_cell_Alexa568_nuc_data(cell_line_id,drug_id8,timepoint_id1,dose_id,signals_647_id6,:,:)),[],1);
    %x2 = reshape((single_cell_Alexa647_nuc_data(cell_line_id,drug_id2,timepoint_id1,dose_id,signals_647_id15,:,:)),[],1);
    %y2 = reshape((single_cell_Alexa568_nuc_data(cell_line_id,drug_id2,timepoint_id1,dose_id,signals_647_id15,:,:)),[],1);
    %x3 = reshape((single_cell_Alexa647_nuc_data(cell_line_id,drug_id5,timepoint_id1,dose_id,signals_647_id15,:,:)),[],1);
    %y3 = reshape((single_cell_Alexa568_nuc_data(cell_line_id,drug_id5,timepoint_id1,dose_id,signals_647_id15,:,:)),[],1);
    %x4 = reshape((single_cell_Alexa647_nuc_data(cell_line_id,drug_id6,timepoint_id1,dose_id,signals_647_id15,:,:)),[],1);
    %y4 = reshape((single_cell_Alexa568_nuc_data(cell_line_id,drug_id6,timepoint_id1,dose_id,signals_647_id15,:,:)),[],1);
    %x5 = reshape((single_cell_Alexa647_nuc_data(cell_line_id,drug_id8,timepoint_id1,dose_id,signals_647_id15,:,:)),[],1);
    %y5 = reshape((single_cell_Alexa568_nuc_data(cell_line_id,drug_id8,timepoint_id1,dose_id,signals_647_id15,:,:)),[],1);

    subplot(1,4,j);
    scatplot(log10(x1),log10(y1),'circles');
    hold on;
    %set(gca,'FontSize',18,'XLim',[2.25 3.5],'YLim',[2.15 2.7],'Box','on');
    colorbar('delete');
    %subplot(5,2,j+2);
    %scatplot(log10(x2),log10(y2),'circles');
    %hold on;
    %set(gca,'FontSize',18,'XLim',[2.35 3.55],'YLim',[2.1 2.8],'Box','on');
    %colorbar('delete');
    %subplot(5,2,j+4);
    %scatplot(log10(x3),log10(y3),'circles');
    %hold on;
    %set(gca,'FontSize',18,'XLim',[2.25 3.5],'YLim',[2.1 2.8],'Box','on');
    %colorbar('delete');
    %subplot(5,2,j+6);
    %scatplot(log10(x4),log10(y4),'circles');
    %hold on;
    %set(gca,'FontSize',18,'XLim',[2.25 2.5],'YLim',[2.1 2.8],'Box','on');
    %colorbar('delete');
    %subplot(5,2,j+8);
    %scatplot(log10(x5),log10(y5),'circles');
    %hold on;
    %set(gca,'FontSize',18,'XLim',[2.2 2.55],'YLim',[2.1 2.9],'Box','on');
    %colorbar('delete');
    
end