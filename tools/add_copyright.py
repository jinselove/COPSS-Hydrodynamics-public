import fileinput
import glob

cfileNames = glob.glob('*.C')

GNU = '// Parallel Finite Element-General Geometry Ewald-like Method.\n// Copyright (C) 2015-2016 Xujun Zhao, Jiyuan Li, Xikai Jiang\n\n// This code is free software; you can redistribute it and/or\n// modify it under the terms of the GNU General Public\n// License as published by the Free Software Foundation; either\n// version 2.1 of the License, or (at your option) any later version.\n\n\n// This code is distributed in the hope that it will be useful,\n// but WITHOUT ANY WARRANTY; without even the implied warranty of\n// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU\n// Lesser General Public License for more details.\n\n\n// You should have received a copy of the GNU General Public\n// License along with this code; if not, write to the Free Software\n// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA\n\n\n'

for cfname in cfileNames:
	with open (cfname,'r') as original : data = original.read()
	with open (cfname,'w') as modified : modified.write(GNU + data)

hfileNames = glob.glob('*.h')

for hfname in hfileNames:
	with open (hfname,'r') as original : data = original.read()
	with open (hfname,'w') as modified : modified.write(GNU + data)

