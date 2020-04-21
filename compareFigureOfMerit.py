#!/usr/bin/env python                                                          #
#                                                                              #
# Autor: Michela Negro, NASA-GSFC.                                             #
# On behalf of the AMEGO Team.                                                 #
#                                                                              #
#------------------------------------------------------------------------------#

"""Compare IRFs from different files
"""

import os
import re
import ast
import argparse
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
sns.set()
sns.set_style("ticks")



__description__ = 'Produce IRFs comparizon plots'


"""Command-line switches.
"""
formatter = argparse.ArgumentDefaultsHelpFormatter
PARSER = argparse.ArgumentParser(description=__description__,
                                 formatter_class=formatter)
PARSER.add_argument('-f', '--infiles', type=str, required=True, nargs='*',
                    help='input .txt file(s) (generated with plotFigureOfMwerit.py)')
PARSER.add_argument('-xs', '--xscale', type=str, choices=['linear', 'log'],
                    default='log', help='specify the scale of the x axis')
PARSER.add_argument('-ys', '--yscale', type=str, choices=['linear', 'log'],
                    default='linear', help='specify the scale of the y axis')
PARSER.add_argument('-ym', '--ymin', type=float, default=None, 
					help='Minimum y-axis value')
PARSER.add_argument('-yM', '--ymax', type=float, default=None, 
					help='Maximum y-axis value')
PARSER.add_argument('-xm', '--xmin', type=float, default=0.1, 
					help='Minimum x-axis value')
PARSER.add_argument('-xM', '--xmax', type=float, default=500,
					help='Maximum x-axis value')
PARSER.add_argument('-t', '--title', type=str, default=' ', 
					help='Plot title')
PARSER.add_argument('-yl', '--ylabel', type=str, default='Figure of Merit', 
					help='labe to display on the y axis')
PARSER.add_argument('-s', '--save', type=ast.literal_eval, choices=[True, False],
                    default=False, help='set to True to save the plot map')
                    
FLAG_STR = '(?<=\_)\w+(?=\.txt)'
LABEL_STR = '^.\w*(?=\_)(?=.*)(?=\_[AE])'
ANG_STR = '(?<=Cos)[0-9]\.[0-9]'

TC_color = np.concatenate((sns.color_palette("Greens_r"), sns.color_palette("Greens_r")))
UC_color = sns.color_palette("Blues_r")
P_color = sns.color_palette("Reds_r")

lines = ['-', '--', '-.', ':', ':', '_.', '--', '-']

def parse_figureofmerit_file(fom_file):
	print('Parsing %s ...'%fom_file)
	file_basename = os.path.basename(fom_file)
	m = re.search(r'%s'%FLAG_STR, file_basename)
	flag = m.group(0)
	m1 = re.search(r'%s'%LABEL_STR, file_basename)
	label = m1.group(0)
	m3 = re.search(r'%s'%ANG_STR, file_basename)
	ang = m3.group(0)
	en, fom = np.loadtxt(fom_file).T
	return en, fom, flag, float(ang), label
	
	
def compare_figureofmerit(**kwargs):
    """plot                                                                                                                                                            
    """
    
    file_list = kwargs['infiles']
    
    TC_count = 0
    UC_count = 0
    P_count = 0
    
    plt.figure(figsize=(10, 7), facecolor='white')
    plt.title(kwargs['title'], size=20)
    for i, f in enumerate(file_list):
        en_, fom_, flag, ang, label = parse_figureofmerit_file(f)
        if flag == 'TC':
        	legend_label = '%s (cos=%.1f) Tracked Compton' %(label, ang)
        	c = TC_color[TC_count]
        	l = lines[TC_count]
        elif flag == 'UC':
        	legend_label = '%s (cos=%.1f) Untracked Compton' %(label, ang)
        	c = UC_color[UC_count]
        	l = lines[UC_count]
        else:
        	legend_label = '%s (cos=%.1f) Pair' %(label, ang)
        	c = P_color[P_count]
        	l = lines[P_count]
        plt.plot(en_, fom_, 'o', color=c, linestyle=l, label=legend_label)
        plt.xlabel('Energy [MeV]', size=18)
        plt.ylabel(kwargs['ylabel'], size=18)
        plt.xlim(kwargs['xmin'], kwargs['xmax'])
        left, right = plt.ylim()
        if kwargs['ymin'] is not None and kwargs['ymax'] is not None:
        	plt.ylim(kwargs['ymin'], kwargs['ymax'])
        elif kwargs['ymin'] is not None and kwargs['ymax'] is None:
        	plt.ylim(kwargs['ymin'], right)
        elif kwargs['ymin'] is None and kwargs['ymax'] is not None:
        	plt.ylim(left, kwargs['ymax'])
        plt.xscale(kwargs['xscale'])
        plt.yscale(kwargs['yscale'])
        plt.legend(fontsize=15)
        plt.grid(False)
        
        if flag == 'TC':
        	TC_count += 1
        elif flag == 'UC':
        	UC_count += 1
        else:
        	P_count += 1
        	
    if kwargs['save']:
    	plt.savefig('Comparison_plot.pdf')
    	print("Created Comparison_plot.pdf ...!")
    	
    plt.show()
	
if __name__ == '__main__':
    args = PARSER.parse_args()
    compare_figureofmerit(**args.__dict__)