% batch function for use with RatCatcher
% collects the mean spike waveforms from a root object
% for a specified filename/filecode

function batchFunction_parallel(bin_id, bin_total, location, batchname, outfile, test)

  if ~test
    addpath(genpath('/projectnb/hasselmogrp/hoyland/RatCatcher'))
    addpath(genpath('/projectnb/hasselmogrp/hoyland/srinivas.gs_mtools'))
    addpath(genpath('/projectnb/hasselmogrp/hoyland/CMBHOME'))
    addpath(genpath('/projectnb/hasselmogrp/hoyland/ln-model-of-mec-neurons'))
  end

  [bin_start, bin_finish] = getParallelOptions(bin_id, bin_total, location, batchname);

  % set up 'local' parallel pool cluster
  pc = parcluster('local');

  % discover the number of available cores assigned by SGE
  nCores = str2num(getenv('NSLOTS'));

  % set up directory for temporary parallel pool files
  parpool_tmpdir = ['~/.matlab/local_cluster_jobs/ratcatcher/ratcatcher_' num2str(bin_id)];
  mkdir(parpool_tmpdir);
  pc.JobStorageLocation = parpool_tmpdir;

  % start parallel pool
  parpool(pc, nCores);

  parfor ii = bin_start:bin_finish

    % set up dummy outfile variable
    outfile_pc = outfile;
    outfile_pc = [outfile_pc(1:end-3), 'mat'];

    % acquire the filename and filecode
    % the filecode should be the "cell number" as a 1x2 vector
    [filename, filecode] = RatCatcher.read(ii, location, batchname);

    % load the data
    % expect a 1x1 Session object named "root"
    % load(filename);
    % root.cel = filecode;

    Results = run_LNmodel_HD_ratCatcher(filename,filecode);

    save(outfile_pc, 'Results');

  end

% end % function
