function export(out,folder)
%EXPORT Portable trajectories and all parameters/options in a MAT file.
if ~exist(folder,'dir'),mkdir(folder);end
% Store callback identity as text: portable MAT v7 cannot encode every handle.
% Restore the documented function handle explicitly before reusing parameters.
for name={'reactionFcn','exchangeFcn'}
 if isa(out.parameters.(name{1}),'function_handle')
  out.parameters.(name{1})=func2str(out.parameters.(name{1}));
 end
end
out.serialization='Callback handles are stored as text; restore before rerunning.';
save(fullfile(folder,'result.mat'),'out','-v7');
n=out.parameters.N;fid=fopen(fullfile(folder,'trajectory.csv'),'w');
assert(fid>=0,'plt:io','Cannot open trajectory output.');cleaner=onCleanup(@()fclose(fid));
fprintf(fid,'time_h,mean_matrix_C,min_matrix_C,max_matrix_C,conditioning_kWh\n');
data=[out.t/3600,mean(out.y(:,n+1:2*n),2),min(out.y(:,n+1:2*n),[],2), ...
 max(out.y(:,n+1:2*n),[],2),out.y(:,end)];
fprintf(fid,'%.12g,%.12g,%.12g,%.12g,%.12g\n',data');
end
