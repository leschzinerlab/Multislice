#!/usr/bin/env python

#Make many simulated stacks for averaging

import subprocess

i = 1 
max = 100

while i <= max:

	cmd = './make_stack_range_defocus.py -i gold_water_new2_256.tif -m ../mtf.spi --min=0 --max=40000 --ds=500 --e=25'
	subprocess.Popen(cmd,shell=True).wait()		

	cmd = 'proc2d image_sim_series_gold_water_new2_256_noise_00000_to_40000.spi image_sim_series_gold_water_new2_256_noise_00000_to_40000invert%03d.spi invert' %(i)
	subprocess.Popen(cmd,shell=True).wait()

	cmd = 'rm image_sim_series_gold_water_new2_256_noise_00000_to_40000.spi'
	subprocess.Popen(cmd,shell=True).wait()

	i = i + 1 
