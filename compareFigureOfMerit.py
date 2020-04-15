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
import matplotlib.pyplot as plt


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
PARSER.add_argument('-xM', '--xmax', type=float, default=100,
					help='Maximum x-axis value')
PARSER.add_argument('-t', '--title', type=str, default=' ', 
					help='Plot title')
PARSER.add_argument('-yl', '--ylabel', type=str, default='Figure of Merit', 
					help='labe to display on the y axis')
PARSER.add_argument('-s', '--save', type=ast.literal_eval, choices=[True, False],
                    default=False, help='set to True to save the plot map')
                    
FLAG_STR = '(?<=\_)\w+(?=\.txt)'
LABEL_STR = '^\S+(?=\_[AE])'

TC_color = 'green'
UC_color = 'blue'
P_color = 'darkred'

def parse_figureofmerit_file(fom_file):
	file_basename = os.path.basename(fom_file)
	m = re.search(r'%s'%FLAG_STR, file_basename)
	flag = m.group(0)
	m1 = re.search(r'%s'%LABEL_STR, file_basename)
	label = m1.group(0)
	en, fom = np.loadtxt(fom_file).T
	return en, fom, flag, label
	
	
def compare_figureofmerit(**kwargs):
    """plot                                                                                                                                                            
    """
    
    file_list = kwargs['infiles']
    
    shade_TC = 1
    shade_UC = 1
    shade_P = 1
    
    plt.figure(figsize=(10, 7))
    plt.title(kwargs['title'], size=20)
    for i, f in enumerate(file_list):
        en_, fom_, flag, label = parse_figureofmerit_file(f)
        if flag == 'TC':
        	legend_label = '%s Tracked Compton' %label
        	c = TC_color
        elif flag == 'UC':
        	legend_label = '%s Untracked Compton' %label
        	c = UC_color
        else:
        	legend_label = '%s Pair' %label
        	c = P_color
        plt.plot(en_, fom_, '-o', color=c, alpha=shade_TC, label=legend_label)
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
        plt.legend(fontsize=18)
        
        if flag == 'TC':
        	shade_TC -= 0.3
        elif flag == 'UC':
        	shade_UC -= 0.3
        else:
        	shade_P -= 0.3
        	
    if kwargs['save']:
    	plt.savefig('Comparison_plot.pdf')
    	
    plt.show()
	
if __name__ == '__main__':
    args = PARSER.parse_args()
    compare_figureofmerit(**args.__dict__)