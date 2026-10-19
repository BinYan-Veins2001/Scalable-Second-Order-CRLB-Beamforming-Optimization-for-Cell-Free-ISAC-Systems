function archive = simtwc_compute_figure(id)
%SIMTWC_COMPUTE_FIGURE Dispatch one compute-only figure experiment.

switch id
    case 1; archive=simtwc_compute_fig01();
    case 2; archive=simtwc_compute_fig02();
    case 3; archive=simtwc_compute_fig03();
    case 4; archive=simtwc_compute_fig04();
    case 5; archive=simtwc_compute_fig05();
    case 6; archive=simtwc_compute_fig06();
    case 7; archive=simtwc_compute_fig07();
    case 8; archive=simtwc_compute_fig08();
    case 9; archive=simtwc_compute_fig09();
    case 10; archive=simtwc_compute_fig10();
    case 11; archive=simtwc_compute_fig11();
    case 12; archive=simtwc_compute_fig12();
    case 13; archive=simtwc_compute_fig13();
    case 14; archive=simtwc_compute_fig14();
    case 15; archive=simtwc_compute_fig15();
    case 16; archive=simtwc_compute_fig16();
    case 17; archive=simtwc_compute_fig17();
    case 18; archive=simtwc_compute_fig18();
    case 19; archive=simtwc_compute_fig19();
    case 20; archive=simtwc_compute_fig20();
    otherwise; error('Supported figure ids are 1:20.');
end
end
