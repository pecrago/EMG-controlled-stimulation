main	run_EMG_control_sim					
						
	calls	load_FDI.m				FDI muscle model parameters
		rampVol_arbStim.m			Executes simulation	
		calls	Contessa_ramp.m			calculates the rates of all SMUs activated in a voluntary excitation ramp 
			stimlevel.m			returns normalized stimulation recruitment level based on prior EMGs
			APbyFrame			calculates all APs in a single frame
			EMGbyFrame			calulates the EMG produced by all APs in a frame
			ForcebyFrame			calculates the force produced by all APs in a frame
			gsfm.m				estimates voluntary emg level in current frame
			signalnoise.m			creates band-limited noise to add to EMG
						
		saves	EMG_controlled_stim_xxx.dat	The entire workspace of rampVol_arbStim
