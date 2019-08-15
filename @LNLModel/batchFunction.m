% batch function for use with RatCatcher
% collects the mean spike waveforms from a root object
% for a specified filename/filecode

% function batchFunction(index, batchname, location, outfile, test)
%
%   if ~test
%     addpath(genpath('/projectnb/hasselmogrp/hoyland/RatCatcher'))
%     addpath(genpath('/projectnb/hasselmogrp/hoyland/srinivas.gs_mtools'))
%     addpath(genpath('/projectnb/hasselmogrp/hoyland/CMBHOME'))
%     addpath(genpath('projectnb/hasselmogrp/hoyland/ln-model-of-mec-neurons'))
%   end
%
%   % acquire the filename and filecode
%   % the filecode should be the "cell number" as a 1x2 vector
%   [filename, filecode] = RatCatcher.read(location, batchname, index);
  filename = 'name';
  filecode = [1, 1];
  % load the data
  % expect a 1x1 Session object named "root"
  load(filename);
  root.cel = filecode;



  %% TODO:
  %% Do whatever needs to be done here to preprocess
  %   you should have the "root" object and the "filecode", so
  %   root.cel = filecode;
  %   and then do whatever
  %% Run the associated functions/scripts
  %% Save the data to a file named "outfile"
  %   e.g. csvwrite(outfile, all_the_data);

  outfile = [outfile(end-2:end) 'mat'];
  save(outfile, 'var1', 'var2', 'var3')

% end % function
