function visualizationDriver()
d = dir('/home/qd4314/Desktop/Figure1and2Update/1stOrder/off/20*');
chooseall = true;

close all


if not(chooseall)
    % choose the last simulation
    folder = d(end-1).name
    prefix = strcat('../out/', folder);
    visualizeOneFolder(prefix)
else % create for all out folders
    for i = 1:length(d)
        folder = d(i).name
        i
        prefix = strcat('/home/qd4314/Desktop/Figure1and2Update/1stOrder/off/', folder);
        try visualizeOneFolder(prefix)
        end
    end
end


save
end

function visualizeOneFolder(prefix)
problemparams = extractParamsFromConfig(prefix);
plottingparams.LW = 3; plottingparams.FS = 40; plottingparams.PP = [0, 0, 18, 18];
% visualizeOverview(prefix, plottingparams)
if problemparams.dim1flag
    visualize1d(prefix, plottingparams, problemparams);
else
    if problemparams.testcaseid == 1 % line source
        visualizeLinesource(prefix, plottingparams, problemparams);
    elseif problemparams.testcaseid == 2 % checkerboard
        visualizeCheckerboard(prefix, plottingparams, problemparams);
    end
end
end
