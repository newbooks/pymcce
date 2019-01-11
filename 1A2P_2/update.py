#!/usr/bin/python

import sys

conformer = {}

fort = sys.argv[1]
headlist = sys.argv[2]

if len(sys.argv) < 2:
	print "Usage: %s fort.38 head_list\n" % sys.argv[0]

for line in open(fort).readlines():
	line = line.strip()
	if line:
		(key, value) = line.split()
		#if value in ['0.000', '1.000']:
		#	conformer[key] = float(value)

		conformer[key] = float(value)

for line in open(headlist).readlines():
	line = line.strip()
	if line:
		tmp = line.split()
		conf_name = tmp[1]
		if conformer.has_key(conf_name):
			new_str = '%s t %3.2f ' % (conf_name, conformer[conf_name])	
			old_str = '%s f %s ' % (conf_name, tmp[3])	
			new_line = line.replace(old_str, new_str)
			print new_line
		else:
			print line

