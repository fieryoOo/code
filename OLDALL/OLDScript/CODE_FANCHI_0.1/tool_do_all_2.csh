#foreach per ( 30 40 60 70 80 90 100 110 120)
foreach per ( 30 40 60 70 80 100 120 )
#csh do_get_good.csh $per
csh tool_do_AtoD_v2_all_1_XP_cv1.csh $per
csh tool_do_surface_slow_lapalce_all.csh $per
csh tool_do_surface_am_laplace.csh $per
csh tool_get_cor_scale.csh $per
csh tool_get_cor_scale_2.csh $per
csh tool_get_scaled_tomography.csh $per
end
