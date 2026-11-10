function study = run_mesh_cost_study(repeats,resume)
%RUN_MESH_COST_STUDY Fixed-tolerance mesh errors and repeated end-to-end CPU costs.
% Wall time includes validation, integration and diagnostics; excludes file I/O.
if nargin<1,repeats=3;end
if nargin<2,resume=false;end
root=fileparts(fileparts(mfilename('fullpath')));names={'heating','phase_cycle','oxygen_limited'};
ns=[10 20 40 80];study=struct('runtime',version,'repeats',repeats,'N',ns,'referenceN',160);
if resume
 old=load(fullfile(root,'results','mesh_cost_report.mat'));study=old.study;
 assert(study.repeats==repeats);
end
% Warm-up avoids counting first function loading only in the smallest run.
[p,s,o]=demo_case('heating',4);o.tEnd=10;plt.simulate(p,s,o);
for j=1:numel(names)
 name=names{j};if resume&&isfield(study,name),continue;end
 [p,s,o]=demo_case(name,160);p.rtol=1e-8;p.atol=1e-11;o.initialStep=1e-6;
 ref=plt.simulate(p,s,o);plt.export(ref,fullfile(root,'results',['reference_' name]));
 times=linspace(0,o.tEnd,25)';rf=interp1(ref.t,ref.y,times);
 rows=zeros(numel(ns),7);timings=zeros(numel(ns),repeats);
 for i=1:numel(ns)
  [p,s,o]=demo_case(name,ns(i));
  for r=1:repeats
   t0=tic;out=plt.simulate(p,s,o);timings(i,r)=toc(t0);
  end
  plt.export(out,fullfile(root,'results',sprintf('mesh_%s_%d',name,p.N)));
  q=interp1(out.t,out.y,times);nb=160/p.N;
  rT=zeros(numel(times),p.N);rX=rT;
  for k=1:p.N
   ii=(k-1)*nb+(1:nb);rT(:,k)=mean(rf(:,160+ii),2);rX(:,k)=mean(rf(:,480+ii),2);
  end
  dT=q(:,p.N+1:2*p.N)-rT;dX=q(:,3*p.N+1:4*p.N)-rX;
  rows(i,:)=[p.N,max(abs(dT(:))),sqrt(mean(dT(:).^2)),max(abs(dX(:))), ...
   median(timings(i,:)),min(timings(i,:)),max(timings(i,:))];
  fprintf('%s N=%d maxT=%.5g C median=%.3f s\n',name,p.N,rows(i,2),rows(i,5));
 end
 item=struct('rows',rows,'timings',timings,'sampleTimes',times,'referenceRuntime',ref.elapsedSeconds);
 study.(name)=item;
 save(fullfile(root,'results','mesh_cost_report.mat'),'study','-v7');
end
end
