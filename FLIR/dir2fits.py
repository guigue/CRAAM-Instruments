#!/usr/bin/python3
import sys, getopt, string, os, glob

def main():

    if (len(sys.argv) < 2) || (len(sys.argv) > 2):
        print(' ')
        print('Usage : ')
        print('   dir2fits <directory path>')
        print(' ')
        print(' ')
        return

    fpfdir = sys.argv[1]

    os.chdir(fpfdir)
    flist = glob.glob('*_fpf')
    for file in flist:
        froot = file[:-4]
        os.spawnl(os.P_WAIT,'/molly/guigue/Programming/CRAAM-Instruments/FLIR/fpf2fits','fpf2zip',froot)
            

if __name__ == "__main__":
    main()
