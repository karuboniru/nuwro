#define PARAMS_ALL()\
PARAM(int,random_seed,0)\
PARAM(int,number_of_events,100000)\
PARAM(int,number_of_test_events,1000000)\
PARAM(int,save_test_events,0)\
PARAM(int,user_events,0)\
PARAM(line,user_params,"")\
PARAM(int,beam_type,0)\
PARAM(line,beam_energy,"1000")\
PARAM(int,beam_particle,14)\
PARAM(vec,beam_direction,"0 0 1")\
PARAM(line,beam_content,"")\
PARAM(string,beam_folder,"flux")\
PARAM(int,beam_file_first,1)\
PARAM(int,beam_file_limit,0)\
PARAM(bool,beam_weighted,0)\
PARAM(vec,beam_offset,"0 0 0")\
PARAM(int,beam_placement,0)\
PARAM(string,beam_inputroot,"")\
PARAM(string,beam_inputroot_flux,"")\
PARAM(string,beam_inputroot_nue,"")\
PARAM(string,beam_inputroot_nueb,"")\
PARAM(string,beam_inputroot_numu,"")\
PARAM(string,beam_inputroot_numub,"")\
PARAM(string,beam_inputroot_nutau,"")\
PARAM(string,beam_inputroot_nutaub,"")\
PARAM(int,beam_test_only,0)\
PARAM(int,target_type,0)\
PARAM(int,nucleus_p,6)\
PARAM(int,nucleus_n,6)\
PARAM(double,nucleus_E_b,34)\
PARAM(double,nucleus_kf,220)\
PARAM(line,target_content,"")\
PARAM(string,geo_file,"target/ND280_v9r7p5.root")\
PARAM(string,geo_name,"ND280Geometry_v9r7p5")\
PARAM(string,geo_volume,"")\
PARAM(vec,geo_o,"0 0 0")\
PARAM(vec,geo_d,"2000 2000 5000")\
PARAM(int,nucleus_target,2)\
PARAM(int,nucleus_model,1)\
PARAM(bool,dyn_qel_cc,1)\
PARAM(bool,dyn_qel_nc,1)\
PARAM(bool,dyn_res_cc,1)\
PARAM(bool,dyn_res_nc,1)\
PARAM(bool,dyn_dis_cc,1)\
PARAM(bool,dyn_dis_nc,1)\
PARAM(bool,dyn_coh_cc,1)\
PARAM(bool,dyn_coh_nc,1)\
PARAM(bool,dyn_mec_cc,0)\
PARAM(bool,dyn_mec_nc,0)\
PARAM(bool,dyn_hip_la,0)\
PARAM(bool,dyn_hip_si,0)\
PARAM(bool,dyn_e_el,0)\
PARAM(bool,dyn_e_spp,0)\
PARAM(double,eel_theta_lab,90)\
PARAM(double,eel_dz,2)\
PARAM(string,eel_alg,"old")\
PARAM(int,qel_kinematics,0)\
PARAM(int,qel_vector_ff_set,2)\
PARAM(int,qel_axial_ff_set,1)\
PARAM(int,qel_rpa,1)\
PARAM(int,qel_strange,1)\
PARAM(int,qel_strangeEM,0)\
PARAM(double,delta_s,0)\
PARAM(double,qel_cc_vector_mass,1000)\
PARAM(double,qel_cc_axial_mass,1030)\
PARAM(double,qel_nc_axial_mass,1030)\
PARAM(double,qel_s_axial_mass,1030)\
PARAM(bool,flux_correction,1)\
PARAM(int,sf_method,0)\
PARAM(bool,sf_fsi,1)\
PARAM(bool,sf_coulomb,1)\
PARAM(int,sf_pb,1)\
PARAM(bool,cc_smoothing,1)\
PARAM(int,delta_FF_set,1)\
PARAM(int,e_spp_ff_set,4)\
PARAM(int,delta_selfenergy,0)\
PARAM(double,pion_axial_mass,0.94)\
PARAM(double,pion_C5A,1.19)\
PARAM(int,delta_angular,2)\
PARAM(int,spp_precision,500)\
PARAM(double,res_dis_cut,1600)\
PARAM(double,bkgrscaling,0.0)\
PARAM(int,res_kind,1)\
PARAM(int,res_hybrid_sampling,1)\
PARAM(int,res_hybrid_resampling,0)\
PARAM(bool,coh_mass_correction,1)\
PARAM(bool,coh_new,1)\
PARAM(int,coh_kind,2)\
PARAM(int,mec_kind,3)\
PARAM(double,mec_ratio_pp,0.85)\
PARAM(double,mec_ratio_ppp,0.8)\
PARAM(double,mec_central_motion,0.0)\
PARAM(double,mec_back_to_back_smearing,0.0)\
PARAM(int,mec_pb_trials,30)\
PARAM(bool,MEC_pauli_blocking,1)\
PARAM(double,MEC_cm_direction,0.0)\
PARAM(bool,kaskada_on,1)\
PARAM(double,kaskada_w,7)\
PARAM(bool,kaskada_redo,0)\
PARAM(bool,kaskada_writeall,0)\
PARAM(string,formation_zone,"fz-new")\
PARAM(double,tau,8)\
PARAM(double,formation_length,1)\
PARAM(bool,first_step,1)\
PARAM(double,step,0.2)\
PARAM(double,kaskada_NN_mfp_scale,1)\
PARAM(int,kaskada_NN_xsec,2)\
PARAM(int,kaskada_NN_inel,2)\
PARAM(int,kaskada_NN_angle,3)\
PARAM(int,kaskada_NN_corr,1)\
PARAM(int,kaskada_piN_xsec,1)\
PARAM(bool,pauli_blocking,1)\
PARAM(bool,mixed_order,1)\
PARAM(double,res_dis_blending_start,2800)\
PARAM(double,res_dis_blending_end,3200)\
PARAM(int,mh_sample_interval,80)\
PARAM(bool,use_weighted_channel,1)\

