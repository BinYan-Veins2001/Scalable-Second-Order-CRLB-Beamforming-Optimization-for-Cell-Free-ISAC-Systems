function test_theorem3_plot()
%TEST_THEOREM3_PLOT Check single-panel CDF and failed-trial accounting.
root=fileparts(fileparts(mfilename('fullpath'))); addpath(genpath(root),'-begin');
s=simtwc_setup();
options=simtwc_theorem3_options(struct('output_tag','plot_test'));
options.regimes=options.regimes(1:2);
trials=cell(2,3);
trials{1,1}=mock_trial('completed',0,0);
trials{1,2}=mock_trial('completed',0,0.5);
trials{1,3}=mock_trial('completed',0,2);
trials{2,1}=mock_trial('completed',0,1);
trials{2,2}=mock_trial('completed',1,0.1);
trials{2,3}=mock_trial('error',0,0.1);
archive=fullfile(s.data_dir,'fig20_plot_test.mat');
save(archive,'options','trials');
fig=simtwc_plot_fig20(archive);
cleanup=onCleanup(@() close(fig)); %#ok<NASGU>
axs=findall(fig,'Type','axes');
assert(numel(axs)==1,'Fig. 20 must have exactly one axes.');
assert(strcmp(axs.XScale,'log'));
assert(strcmp(axs.XLabel.Interpreter,'latex'));
assert(strcmp(axs.YLabel.Interpreter,'latex'));
assert(abs(fig.Position(3)-3.5)<1e-12);
u=fig.UserData;
assert(isequal(u.joint_ratio(1,:),[0 0.5 2]));
assert(isequal(u.joint_ratio(2,:),[1 inf inf]));
assert(norm(u.attainment_rate-[2/3;1/3],inf)<1e-12);
curves=findall(axs,'Type','stair');
assert(numel(curves)==2);
for j=1:numel(curves)
    assert(all(diff(curves(j).YData)>=0));
    assert(all(curves(j).XData>0));
end
fprintf('THEOREM 3 SINGLE-PANEL PLOT TEST PASSED.\n');
end

function trial=mock_trial(status,primal,ratio)
trial=struct('status',status,'audit',struct('primal',primal,'joint_ratio',ratio), ...
    'solution',struct('terminal_options',struct('feas_tol',1e-7)));
end
