# encoding: utf-8
'''
A tool for generating model configuration files...

Created on 8 Feb 2019

@author: Eugen Rožić
'''
import sys, os, re


template_path = sys.argv[1]

out_folder = template_path+'_out'
if not os.path.exists(out_folder):
    os.makedirs(out_folder)

version = int(template_path[-1])

SSs = [0.0] #3.25
SB_tips = [1.0]
SB_crosss = [4.5] #absolute value, SB_tip/2 will be deduced
SB_centrals = [1.5]
BB_centrals = [6.0] #relative to SB_centrals
SB_sides = [1.5]
BB_sides = [2.0, 4.0] #relative to SB_sides

D_mu = 15.0
VX = 5.0

def calc_central_int(version, BB_central, SB_central, SB_cross, SB_tip):
    if version == 2:
        return 6.0*BB_central - 7.0*SB_central - 2*SB_cross - SB_tip
    elif version == 3:
        return 6.5*BB_central - 7.0*SB_central - 2*SB_cross - SB_tip
    else:
        raise Exception()
    
def calc_side_int(version, BB_side, SB_side, SB_cross, SB_tip):
    if version == 2 or version == 3:
        return 6*(BB_side-SB_side) - 0.4*SB_cross - SB_tip
    else:
        raise Exception()

def match_replace(matchobj):
    return str(globals()[matchobj.group(1)])

def fill_template(line):
    return re.sub(r'<([a-zA-Z_][a-zA-Z_0-9\-]*)>', match_replace, line)

for SS in SSs:
    if SS == 0:
        if_SS_zero = "(VX, vx)#"
    else:
        if_SS_zero = ""
    for SB_tip in SB_tips:
        if SB_tip == 0:
            if_SB_tip_zero = "(VX, vx)#"
        else:
            if_SB_tip_zero = ""
        for SB_cross in SB_crosss:
            SB_cross -= SB_tip/2
            for SB_central in SB_centrals:
                for BB_central in BB_centrals:
                    BB_central += SB_central
                    for SB_side in SB_sides:
                        for BB_side in BB_sides:
                            BB_side += SB_side
                            
                            cent_eff = calc_central_int(version, BB_central, SB_central, SB_cross, SB_tip)
                            side_eff = calc_side_int(version, BB_side, SB_side, SB_cross, SB_tip)
                            
                            cfg_filename = "{:.2f}-{:.2f}-{:.2f}-{:.2f}_{:.2f}-{:.2f}_{:.2f}.cfg".format(
                                SS, SB_tip, SB_cross, SB_central, BB_central, SB_side, BB_side)
                            cfg_filepath = os.path.join(out_folder, cfg_filename)
                            
                            with open(cfg_filepath, 'w') as cfg_file:
                                with open(template_path, 'r') as template_file:
                                    for line in template_file:
                                        cfg_file.write(fill_template(line))
                            
