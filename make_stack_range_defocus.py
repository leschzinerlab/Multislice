#!/usr/bin/env python


import optparse
from sys import *
import os,sys,re
from optparse import OptionParser
import glob
import subprocess
from os import system
import linecache

def setupParserOptions():
        parser = optparse.OptionParser()
        parser.set_usage("%prog -i <image.tif> -m <mtf> --min=<float> --max=<float> --ds=<float> --e=<float>")
        parser.add_option("-i",dest="image",type="string",metavar="FILE",
                help=".tif image from multislice)")
        parser.add_option("-m",dest="mtf",type="string",metavar="FILE",
                help="SPIDER text file with MTF-info for detector")
        parser.add_option("--min",dest="min",type="int", metavar="INT",
                help="Minimum defocus for range (1 um = 10,000 A)")
        parser.add_option("--max",dest="max",type="int", metavar="INT",
                help="Maximum defocus for range (1 um = 10,000 A)")
        parser.add_option("--ds",dest="step",type="int", metavar="INT",
                help="Step size for range of defocuses")
        parser.add_option("--e",dest="dose",type="int", metavar="INT",
                help="Electron dose (e/A^2)")	
        parser.add_option("-n", action="store_true",dest="n",default=False,
                help="Flag for no noise added to images")
        parser.add_option("-d", action="store_true",dest="debug",default=False,
                help="debug")
        options,args = parser.parse_args()

        if len(args) > 0:
                parser.error("Unknown commandline options: " +str(args))

        if len(sys.argv) < 4:
                parser.print_help()
                sys.exit()
        params={}
        for i in parser.option_list:
                if isinstance(i.dest,str):
                        params[i.dest] = getattr(options,i.dest)
        return params

def main(params):

	i = params['image']
	mtf = params['mtf']
	MIN = int(params['min'])
	MAX = int(params['max'])
	step = int(params['step'])
	dose = params['dose']
	n = params['n']
	debug = params['debug']
	out = '%s_defocus.par' %(i[:-4])
	o1 = open(out,'w')

	l = MIN

	while l <= MAX:
	
		if n is True:
			cmd = '/labdata/allab/michaelc/Scripts/Multislice/generate_gold_images.py -i %s -m %s --def=%s --e=%s -n' %(i,mtf,str(l),dose)
			subprocess.Popen(cmd,shell=True).wait()
	
			i2 = '%s_%05dA_%02ddose_image.spi'%(i[:-4],(l),dose)
			if debug is True:
				print 'proc2d %s image_sim_series_%s.spi' %(i2,i[:-4])

			o1.write('%s\n' %(l))			

			cmd = 'proc2d %s image_sim_series_%s_%05d_to_%05d.spi' %(i2,i[:-4],MIN,MAX)
			subprocess.Popen(cmd,shell=True).wait()

			cmd = 'rm %s_*' %(i[:-4])
			subprocess.Popen(cmd,shell=True).wait()

		if n is False:
			cmd = '/labdata/allab/michaelc/Scripts/Multislice/generate_gold_images.py -i %s -m %s --def=%s --e=%s' %(i,mtf,str(l),dose)
			subprocess.Popen(cmd,shell=True).wait()

			i2 = '%s_noise_%05dA_%02ddose_image.spi'%(i[:-4],(l),dose)

			o1.write('%s\n' %(l))

                        cmd = 'proc2d %s image_sim_series_%s_noise_%05d_to_%05d.spi' %(i2,i[:-4],MIN,MAX)
			subprocess.Popen(cmd,shell=True).wait()

                        cmd = 'rm %s_*' %(i[:-4])
                        subprocess.Popen(cmd,shell=True).wait()
				

		l = l + step
	

if __name__ == "__main__":     
        params=setupParserOptions()     
        main(params)
