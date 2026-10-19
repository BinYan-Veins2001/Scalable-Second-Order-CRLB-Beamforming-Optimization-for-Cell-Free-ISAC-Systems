function [geo,g,c,G,out]=paperfig_geometry(p,z,nu,options)
%PAPERFIG_GEOMETRY Re-evaluate the current paper's coordinated NNLS geometry.
if nargin<4; options=struct(); end
o=costar_paper_options(options); N=p.cfg.sys.M*p.cfg.sys.Nt; U=p.cfg.sys.U;
W=real_vector_to_complex_matrix(z,N,U);
[~,gc,out]=compute_crlb_objective_grad(p.cfg,p.smat,W);
g=complex_matrix_to_real_vector(gc);
[c,G]=costar_sosp_constraints(p.cfg,p.channels,W);
geo=costar_paper_geometry(g,c,G,U,nu,o);
end
