# EQ36parser
Python script to clean output from EQ6 model runs

This python script takes a tab file created by the EQ6 model, parses the file, and creates an excel spreadsheet with information important to most users. In addition, it saves three plots:
  modes.png - plots the molar modal abundance of mineral phases versus logzi
  minerals.png - plots the logmoles of minerals versus logzi
  Fe23.png - plots the Fe3+/Fetot ratio in the solid rock (computed assuming all Fe3+ in phase magnetite) versus logzi

CURRENT SUPPORTED MINERAL PHASES:
	MAGNETITE
	BRUCITE
	CLINOCHLORE
	TREMOLITE
	WOLLASTONITE
	OLIVINE (solid solution)
	ORTHOPYROXENE (solid solution)
	CLINOPYROXENE (solid solution)
	BIOTITE (solid solution)
	GARNET (solid solution)

Known issues:
  Not all mineral phases are inlcuded. The script may error if a mineral phase is produced that is not recognized by the script. If you     are using this script and would like a mineral phase added, please contact Kayla Iacovino at kayla.iacovino@asu.edu.
