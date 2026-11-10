function startup_porousltne()
%STARTUP_POROUSLTNE Add source, examples and tests without recursive paths.
root=fileparts(mfilename('fullpath'));
addpath(fullfile(root,'src'),fullfile(root,'examples'),fullfile(root,'tests'));
fprintf('PorousLTNE 0.2.0 initialized. Run example_compost or test_porousltne.\n');
end
