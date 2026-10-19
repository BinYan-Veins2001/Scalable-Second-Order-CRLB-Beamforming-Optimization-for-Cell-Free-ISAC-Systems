function fig = simtwc_plot_figure(id)
%SIMTWC_PLOT_FIGURE Dispatch one disk-only plot and export routine.

switch id
    case 1; fig=simtwc_plot_fig01();
    case 2; fig=simtwc_plot_fig02();
    case 3; fig=simtwc_plot_fig03();
    case 4; fig=simtwc_plot_fig04();
    case 5; fig=simtwc_plot_fig05();
    case 6; fig=simtwc_plot_fig06();
    case 7; fig=simtwc_plot_fig07();
    case 8; fig=simtwc_plot_fig08();
    case 9; fig=simtwc_plot_fig09();
    case 10; fig=simtwc_plot_fig10();
    case 11; fig=simtwc_plot_fig11();
    case 12; fig=simtwc_plot_fig12();
    case 13; fig=simtwc_plot_fig13();
    case 14; fig=simtwc_plot_fig14();
    case 15; fig=simtwc_plot_fig15();
    case 16; fig=simtwc_plot_fig16();
    case 17; fig=simtwc_plot_fig17();
    case 18; fig=simtwc_plot_fig18();
    case 19; fig=simtwc_plot_fig19();
    case 20; fig=simtwc_plot_fig20();
    otherwise; error('Supported figure ids are 1:20.');
end
end
