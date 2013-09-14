foreach net (TA XR)
foreach period (50 80)
foreach event (20080830065405 20081028160003 20081002232809 20091011031213 20100113050257)
Plot_one_event.csh /home/tianye/data_Eikonal/SAC_$net $period $event /home/tianye/data_Eikonal/SAC_XR/region_XR
end
end
end
