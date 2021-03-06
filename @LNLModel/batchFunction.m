% batch function for use with RatCatcher
% collects the mean spike waveforms from a root object
% for a specified filename/filecode

function batchFunction(index, location, batchname, outfile, test)

  if ~test
    addpath(genpath('/projectnb/hasselmogrp/hoyland/RatCatcher'))
    addpath(genpath('/projectnb/hasselmogrp/hoyland/srinivas.gs_mtools'))
    addpath(genpath('/projectnb/hasselmogrp/hoyland/CMBHOME'))
    addpath(genpath('/projectnb/hasselmogrp/hoyland/ln-model-of-mec-neurons'))
  end

  % acquire the filename and filecode
  % the filecode should be the "cell number" as a 1x2 vector
  [filename, filecode] = RatCatcher.read(index, location, batchname);

  % load the data
  % expect a 1x1 Session object named "root"
  % load(filename);
  % root.cel = filecode;

  Results = run_LNmodel_HD_ratCatcher(filename,filecode);

  outfile = [outfile(1:end-3) 'mat'];
  save(outfile, 'Results');

% end % function
