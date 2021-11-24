function convergenceDriver()
d = dir('/home/thomas/Work/Projects/rSN/code/out/save_2d_lowrank_vs_fullrank_CB/run2/20*');
chooseall = true;
close all

rhosCB = {};
if not(chooseall)
    % choose the last simulation
    folder = d(end-1).name
    prefix = strcat('../out/', folder);
    rhosCB{end+1} = visualizeOneFolder(prefix)
else % create for all out folders
    for i = 1:length(d)
        folder = d(i).name
        prefix = strcat('/home/thomas/Work/Projects/rSN/code/out/save_2d_lowrank_vs_fullrank_CB/run2/', folder);
        rhosCB{end+1} = visualizeOneFolder(prefix)
    end
end


save('rhosCB','rhosCB')

end

function rho = visualizeOneFolder(prefix)
problemparams = extractParamsFromConfig(prefix);
plottingparams.LW = 3; plottingparams.FS = 40; plottingparams.PP = [0, 0, 18, 18];
visualizeOverview(prefix, plottingparams)
if problemparams.dim1flag
    rho = visualize1d(prefix, plottingparams, problemparams);
else
    if problemparams.testcaseid == 1 % line source
        rho = visualizeLinesource(prefix, plottingparams, problemparams);
    elseif problemparams.testcaseid == 2 % checkerboard
        rho = visualizeCheckerboard(prefix, plottingparams, problemparams);
    end
end
end
