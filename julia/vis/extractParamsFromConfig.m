function params = extractParamsFromConfig(folder)

filename = strcat(folder, '/config.txt');
f = fileread(filename);

nx = regexp(f, '(?<=nx = )\w*', 'match'); % matches the word after 'nx = '
params.nx = str2num(nx{1});

ny = regexp(f, '(?<=ny = )\w*', 'match');
params.ny = str2num(ny{1});

cfl = regexp(f, '(?<=cfl = )\w*', 'match');
params.cfl = str2num(cfl{1});

quadraturetype = regexp(f, '(?<=quadraturetype = )\w*', 'match');
params.quadraturetype = str2num(quadraturetype{1});

quadratureorder = regexp(f, '(?<=quadratureorder = )\w*', 'match');
params.quadratureorder = str2num(quadratureorder{1});


convolutionStrength = regexp(f, '(?<=convolutionmagnitude = )\w*', 'match');
params.convolutionStrength = str2num(convolutionStrength{1});
if params.convolutionStrength==0
    params.convolutionflag = false;
else
    params.convolutionflag = true; 
end

rotationmagnitude = regexp(f, '(?<=rotationmagnitude = )\w*', 'match');
params.rotationmagnitude = str2num(rotationmagnitude{1});

testcaseid = regexp(f, '(?<=testcaseid = )\w*', 'match');
params.testcaseid = str2num(testcaseid{1});

dim1flag = regexp(f, '(?<=dim1flag = )\w*', 'match');
params.dim1flag = str2num(dim1flag{1});

whichrank = regexp(f, '(?<=whichrank = )\w*', 'match');
params.whichrank = str2num(whichrank{1});
if params.whichrank == 0
    params.lowrankflag = false;
else
    params.lowrankflag = true;
end

if params.quadraturetype == 1
    nq = params.quadratureorder^2;
elseif params.quadraturetype == 2
    nq = 4*params.quadratureorder^2-8*params.quadratureorder+6;
elseif params.quadraturetype == 3
    nq = 10*params.quadratureorder^2-20*params.quadratureorder+12;
elseif params.quadraturetype == 4
    nq = params.quadratureorder;
else
    warning('What quadrature type is this?');
    nq = 0;
end
params.nquadpoints = nq;

end
