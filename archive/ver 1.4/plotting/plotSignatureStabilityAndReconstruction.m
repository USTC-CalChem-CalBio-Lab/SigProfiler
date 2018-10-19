function plotSignatureStabilityAndReconstruction(x, stability, reconstruction, input)

% Ludmil B. Alexandrov
% Cancer Genome Project
% Wellcome Trust Sanger Institute
% la2@sanger.ac.uk
%
% This software and its documentation are copyright 2012 by the
% Wellcome Trust Sanger Institute/Genome Research Limited. All rights are reserved.
% This software is supplied without any warranty or guaranteed support whatsoever. 
% Neither the Wellcome Trust Sanger Institute nor Genome Research Limited 
% is responsible for its use, misuse, or functionality.

    screen_size = get(0, 'ScreenSize');
    f1 = figure('InvertHardcopy','off','Color',[1 1 1]);
    set(f1, 'Position', [0 0 screen_size(3) screen_size(4) ] );
    [AX, H1, H2] = plotyy(x, stability, x, reconstruction, 'plot');
    %set(AX(1), 'YGrid', 'on');
    set(AX(1), 'YTick',[0.4 0.5 0.6 0.7 0.8 0.85 0.9 0.95 1.0]);
    set(AX(1), 'ylim',[0.4 1.01]);

    set(get(AX(1),'Ylabel'),'String','Signature Stability', 'FontSize', 30);
    set(get(AX(2),'Ylabel'),'String','Average Frobenius Reconstruction Error', 'FontSize', 30) 
    set(AX,{'FontSize'},{20;20});  
    set(AX,{'FontWeight'},{'bold';'bold'});
    set(AX,{'LineWidth'},{3;3});
    set(AX,{'ycolor'},{[0.847058832645416 0.160784319043159 0];[0.0392156876623631 0.141176477074623 0.415686279535294]});
    set(H1, 'MarkerFaceColor',[1 0.600000023841858 0.7843137383461],...
        'MarkerEdgeColor',[0.847058832645416 0.160784319043159 0],...
        'MarkerSize',20,...
        'Marker','o',...
        'LineWidth',2,...
        'LineStyle',':',...
        'Color',[0.847058832645416 0.160784319043159 0]);
    set(H2, 'MarkerFaceColor',[0.0392156876623631 0.141176477074623 0.415686279535294],...
        'MarkerEdgeColor',[0.701960802078247 0.780392169952393 1],...
        'MarkerSize',20,...
        'Marker','square',...
        'LineWidth',2,...
        'LineStyle',':',...
        'Color',[0.0392156876623631 0.141176477074623 0.415686279535294]);

    xlabel('Number of Signatures', 'FontWeight', 'bold', 'FontSize', 30); 
    title(['Signatures in ' num2str(size(input.originalGenomes, 2)) input.cancerType ' samples'],'FontWeight','bold','FontSize',30);
end