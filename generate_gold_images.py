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
        parser.set_usage("%prog -i <image.tif> -m <mtf> --def=<float> --e=<float>")
        parser.add_option("-i",dest="image",type="string",metavar="FILE",
                help=".tif image from multislice)")
        parser.add_option("-m",dest="mtf",type="string",metavar="FILE",
                help="SPIDER text file with MTF-info for detector")
        parser.add_option("--def",dest="def",type="int", metavar="INT",
                help="Defocus (in Angstroms) for output image (1 um = 10,000 A)")
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

def applyMTF(i,mtf):
        spifile = "currentSpiderScript.spi"
        if os.path.isfile(spifile):
                os.remove(spifile)
        spi=open(spifile,'w')
        spicmd="FT\n"
        spicmd+="%s\n" %(i)
        spicmd+="%s_ft\n" %(i)
        spicmd+="FD\n" 
        spicmd+="%s_ft\n" %(i)
        spicmd+="%s_ft_mtf\n" %(i)
        spicmd+="%s\n" %(mtf)
        spicmd+="FT\n" 
        spicmd+="%s_ft_mtf\n" %(i)
        spicmd+="%s_ft_mtf_ft\n" %(i)
        runSpider(spicmd)

def runSpider(lines):
       spifile = "currentSpiderScript.spi"
       if os.path.isfile(spifile):
               os.remove(spifile)
       spi=open(spifile,'w')
       spi.write("MD\n")
       spi.write("TR OFF\n")
       spi.write("MD\n")
       spi.write("VB OFF\n")
       spi.write("MD\n")
       spi.write("SET MP\n")
       spi.write("(0)\n")
       spi.write("\n")
       spi.write(lines)

       spi.write("\nEN D\n")
       spi.close()
       spicmd = "spider spi @currentSpiderScript"
       spiout = subprocess.Popen(spicmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).stderr.read()
       output = spiout.strip().split()
       if "ERROR" in output:
               print "Spider Error, check 'currentSpiderScript.spi'\n"
               sys.exit()
       # clean up
       os.remove(spifile)
       if os.path.isfile("LOG.spi"):
               os.remove("LOG.spi")
       resultf = glob.glob("results.spi.*")
       if resultf:
               for f in resultf:
                       os.remove(f)

def main(params):
	debug = params['debug']	
	image = params['image']
	dose = params['dose']
	defocus = int(params['def'])
	mtf = params['mtf']
	n = params['n']
	if n is False:

		#Add noise to .tif image:  add_poisson
		cmd="/labdata/allab/michaelc/Scripts/Multislice/bin/add_poisson"
		
		if debug is True:
			print "%s" %(cmd)
			print "%s" %(image)
			print "%s_noise.tif" %(image[:-4])
			print "%s" %(dose)

		a = subprocess.Popen(cmd, stdin=subprocess.PIPE,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		[o,e] = a.communicate("%s\n"%(image) + "%s_noise.tif\n"%(image[:-4]) + "%s\n"%(dose))
		i3 = "%s_noise.tif"%(image[:-4])

	else:
		i3 = "%s"%(image)

	#Simulate EM image: image
	
	cmd2="/labdata/allab/michaelc/Scripts/Multislice/bin/image"
	if debug is True:
		print "%s" %(cmd2)
		print "%s" %(i3)
		print "0"
		print "%s_%05dA.tif" %(i3[:-4],defocus)
		print "2.2"
		print "%s" %(defocus)
		print "10"
		print "0"
		print "0"
		print "0"
		print "0"
		print "0"
		print "0"

        a = subprocess.Popen(cmd2, stdin=subprocess.PIPE,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        [o,e] = a.communicate("%s\n"%(i3) + "0\n" + "%s_%05dA.tif\n"%(i3[:-4],defocus) + "2.2\n" + "%s\n"%(defocus) + "10\n" + "0\n" + "0\n" + "0\n" + "0\n" + "0\n" + "0\n")	

	if debug is True:
		print "/labdata/allab/michaelc/Scripts/Multislice/tif2spi.b %s_%05dA.tif"%(i3[:-4],defocus)

	#Convert .tif to .spi
	cmd = "/labdata/allab/michaelc/Scripts/Multislice/tif2spi.b %s_%05dA.tif"%(i3[:-4],defocus)
	subprocess.Popen(cmd,shell=True).wait()	
	
	#Run spider script:
	i2 = "%s_%05dA"%(i3[:-4],defocus)
	if debug is True:
		print "applyMTF(%s,%s)"%(i2,mtf[:-4])

	applyMTF(i2,mtf[:-4])

	#Normalize image with proc2d
	cmd = "rm %s_%05dA_ft.spi %s_%05dA_ft_mtf.spi" %(i3[:-4],defocus,i3[:-4],defocus)
	subprocess.Popen(cmd,shell=True).wait()

	cmd = "proc2d %s_%05dA_ft_mtf_ft.spi %s_%05dA_%02ddose_image.spi" %(i3[:-4],defocus,i3[:-4],defocus,dose)
	subprocess.Popen(cmd,shell=True).wait()

	cmd = "rm %s_%05dA_ft_mtf_ft.spi" %(i3[:-4],defocus)
	subprocess.Popen(cmd,shell=True).wait()

if __name__ == "__main__":     
        params=setupParserOptions()     
        main(params)
